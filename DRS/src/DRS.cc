#include "DRS.h"
#include <string>
#include <vector>

#include <iostream>
#include <iomanip>
#include <atomic>
void parse_input_argv_DRS(Parameter & param, int argc, char * argv[], 
	std::string & A_filename, std::string & B_filename, std::string & C_filename, Data & data)
{
	for (int i = 1; i < argc; i++)
	{
		if (argv[i][0] != '-')
		{
			break;
		}
		if (++i >= argc)
		{
			//error
		}
		else if (std::string(argv[i - 1]) == "-epoch")
		{
			param.max_iter = atoi(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-nthread")
		{
			param.total_num_threads = atoi(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-block_size")
		{
			param.block_size = atoi(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-check_step")
		{
			param.check_step = atoi(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-eta")
		{
			param.eta = atof(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-lambda")
		{
			param.lambda = atof(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-scale")
		{
			param.scale = atoi(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-choice")
		{
			param.choice = atoi(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-A")
		{
			A_filename = std::string(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-b")
		{
			B_filename = std::string(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-c")
		{
			C_filename = std::string(argv[i]);
		}
		else
		{
			// error
		}
	}
	return;
}

void pre_computing_DRS(Data & data, Parameter & param)
{
	int A_col = param.A_col;
	int A_row = param.A_row;
	data.M.resize(A_row, A_col);
	data.d.resize(A_row);
	
	if (param.scale)
	{
		double scale;

		// left scale
		for (int i = 0; i < A_row; i++)
		{
			scale = data.A.row(i).dot(data.A.row(i));
			scale = 1. / (scale > 0 ? scale : 1);
			scale = sqrt(scale);
			data.A.row(i) *= scale;
			data.b(i) *= scale;
		}
	}
	data.c_linf = data.c.lpNorm<Eigen::Infinity>();
	data.b_linf = data.b.lpNorm<Eigen::Infinity>();
	// block info
	data.block.resize(0);
	for (int block_start = 0; block_start < A_col; block_start+=param.block_size)
	{
		data.block.push_back(block_start);
	}
	param.block_num = data.block.size();
	data.block.push_back(A_col);
	
	// pre-compute
	MatrixXd T(A_row, A_row);
	//data.T = (data.A * data.A.transpose()).inverse();
	T = (data.A * data.A.transpose()).inverse();
	data.T = T;

	//T = T.inverse();
	
	//data.M = T;
	
	//std::cout << (data.M-data.T).trace() << std::endl;
	//char res;
	//std::cin >> res;

	MatrixXd temp(A_col, A_row);
	temp = data.A.transpose() * T;
	data.M = temp * data.A;
	data.d = temp * data.b;
}

void print_parameters(Parameter & param)
{
	using namespace std;
	cout << "Parameter settings:" << endl;
	cout << "---------------------------------" << endl;
	cout.width(28);
	cout << left << "Problem size:" << right << setw(3) << param.A_col << endl;
	cout.width(28);
	cout << left << "eta:" << right << setw(3) << param.eta << endl;
	cout << left << "lambda:" << right << setw(3) << param.lambda << endl;
	cout << "---------------------------------" << endl;
}
using namespace std;

DRS::DRS(VectorXd * z_, VectorXd * s_, VectorXd *Ax_, VectorXd *ATs_, VectorXd *mw_, Data * data_, Parameter * param_,
	StoppingCriteria * crit_, Result * res_)
{
	z = z_;
	s = s_;
	Ax = Ax_;
	ATs = ATs_;
	mw = mw_;
	data = data_;
	param = param_;
	crit = crit_;
	res = res_;
	A_col = param->A_col;
	A_row = param->A_row;
	x.resize(A_col);
	w.resize(A_col);
	u.resize(A_col);
	block_num = param->block_num;
}

DRS::~DRS()
{
}

extern std::atomic<int> block_id;
int DRS::get_id()
{
	return (block_id++) % (block_num);
}

void DRS::update(int id)
{
	x = (*z);
	w = x - param->lambda*data->c;
	w = 0.5*w+0.5*w.cwiseProduct(w.cwiseSign());
	int size = data->block[id + 1] - data->block[id];
	VectorXd delta;
	delta = data->M.block(data->block[id], 0, size, A_col) * (x - 2 * w) - x.segment(data->block[id], size) + w.segment(data->block[id], size) + data->d.segment(data->block[id], size);
	// update z
	(*z).segment(data->block[id],size)+= delta*param->eta;
}

void DRS::parallel_getmin(int id)
{
	int worker_num = param->total_num_threads;
	int size = id < worker_num - 1 ? A_col / worker_num : A_col - (worker_num - 1)*(A_col / worker_num);
	int start = (A_col / worker_num)*id;
	VectorXd temp=(*z).segment(start,size) - param->lambda*data->c.segment(start,size);
	(*mw).segment(start,size) = 0.5*temp+0.5*temp.cwiseProduct(temp.cwiseSign());
}
void DRS::parallel_stop(int id)
{
	int worker_num = param->total_num_threads;
	// copy z
	// norm(Ax-b)
	int size = id < worker_num - 1 ? A_row / worker_num : A_row - (worker_num - 1)*(A_row / worker_num);
	int start = (A_row / worker_num)*id;
	// Ax.resize(size);
	(*Ax).segment(start, size) = data->A.block(start, 0, size, A_col)*(*mw) - data->b.segment(start, size);
	(*ATs).segment(start, size) = data->M.block(start, 0, size, A_col)*(*z) - data->d.segment(start, size) + (*mw).segment(start, size) - (*z).segment(start, size);
	// update s
	(*s).segment(start, size) << data->T.block(start, 0, size, A_row)*(data->A*(*z) - data->b);
	start = (A_row / worker_num)*id;
	size = id < worker_num - 1 ? A_col / worker_num : A_col - (worker_num - 1)*(A_col / worker_num);
}

void DRS::check_stop()
{
	VectorXd temp;
	double primal = (*mw).dot(data->c);
	double dual = (*s).dot(data->b)/param->lambda;
	crit->pfeas = (*Ax).lpNorm<Eigen::Infinity>();
	crit->dfeas = (*ATs).lpNorm<Eigen::Infinity>();
	crit->relgap = fabs(primal - dual) / (1 + fabs(primal) + fabs(dual));
	if (crit->pfeas < param->feastol * (1 + data->b_linf) && crit->dfeas < param->feastol * (1+ data->c_linf) && crit->relgap < param->feastol) 
	{
		res->status = 1;
		param->stop = 1;
		res->optimal = primal;
	}
	else
	{
		crit->pfeas = 0;
		crit->dfeas = 0;
	}
}

void DRS::get_minizer()
{
	for (int i = 0; i < A_col; i++)
	{
		(*z)(i) -= param->lambda*data->c(i);
		(*z)(i) = (*z)(i) > 0 ? (*z)(i) : 0;
	}
}
