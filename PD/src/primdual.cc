#include "primdual.h"
#include <string>
#include <vector>
#include <ctime>
#include <mutex>
#include <iostream>
#include <iomanip>
#include <atomic>
#include <fstream>
#include <string>
using namespace std;
using namespace Eigen;
void parse_input_argv(Parameter & param, int argc, char * argv[], std::string & A_filename, std::string & B_filename, std::string & C_filename,  Data & data)
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
		else if (std::string(argv[i - 1]) == "-block_size_p")
		{
			param.block_size_p = atoi(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-block_size_d")
		{
			param.block_size_d = atoi(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-check_step")
		{
			param.check_step = atoi(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-delta")
		{
			param.delta = atof(argv[i]);
		}
		else if (std::string(argv[i - 1]) == "-alpha")
		{
			param.alpha = atof(argv[i]);
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

void pre_computing_PD(Data & data, Parameter & param)
{
	int A_col = param.A_col;
	int A_row = param.A_row;

	if (param.scale)
	{
		// scale primitive data
		double MinScale = 1e-3;
		double MaxScale = 1e3;
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
	data.b_linf = data.b.lpNorm<Eigen::Infinity>();
	data.c_linf = data.c.lpNorm<Eigen::Infinity>();
	// block info
	data.block.resize(0);
	int block_start;
	for (block_start = 0; block_start < A_row; block_start += param.block_size_d)
	{
		data.block.push_back(block_start);
	}
	if (block_start > A_row)
	{
		block_start = A_row;
	}
	for (; block_start < A_row + A_col; block_start += param.block_size_p)
	{
		data.block.push_back(block_start);
	}
	param.block_num = data.block.size();
	data.block.push_back(A_col+A_row);

	//pre-conditioning
	data.gamma = VectorXd::Zero(param.A_row);
	data.eta = VectorXd::Zero(param.A_col);
	for (int j = 0; j < param.A_col; j++)
	{
		for (int i = 0; i < param.A_row; i++)
		{
			data.gamma(i) += pow(fabs(data.A(i, j)), param.alpha);
			data.eta(j) += pow(fabs(data.A(i, j)), 2.0 - param.alpha);
		}
	}
	for (int j = 0; j < param.A_col; j++)
	{
		data.eta(j) = 1.0 / data.eta(j);
	}
	for (int i = 0; i < param.A_row; i++)
	{
		data.gamma(i) = 1.0 / data.gamma(i);
	}
}

void print_parameters(Parameter & param)
{
	using namespace std;
	cout << "Parameter settings:" << endl;
	cout << "---------------------------------" << endl;
	cout.width(28);
	cout << left << "Problem size:" << right << setw(3) << param.A_col << endl;
	cout.width(28);
	cout << left << "alpha:" << right << setw(3) << param.alpha << endl;
	cout << left << "delta:" << right << setw(3) << param.delta << endl;
	cout << "---------------------------------" << endl;
}

PD::PD(VectorXd * s_, VectorXd * x_, VectorXd * Ax_, VectorXd * ATs_, Data * data_, Parameter * param_,
	StoppingCriteria * crit_, Result * res_)
{
	s = s_;
	x = x_;
	Ax = Ax_;
	ATs = ATs_;
	data = data_;
	param = param_;
	crit = crit_;
	res = res_;
	A_col = param->A_col;
	A_row = param->A_row;
	block_num = param->block_num;
}
extern std::atomic<int> block_id;
int PD::get_id()
{
	return (block_id++) % (block_num);
	//return rand() % block_num;
}
extern std::mutex mu;
void PD::update(int id, int epoch)
{
	double adapt = 1.0/sqrt(epoch);
	//double adapt = 1.0;
	if (data->block[id] < param->A_row)
	{
		VectorXd T;
		int start = data->block[id];
		int size = data->block[id + 1] - data->block[id];
		T = data->gamma.segment(start, size).cwiseProduct((*Ax).segment(start, size) - data->b.segment(start, size));
		(*s).segment(data->block[id], data->block[id + 1] - data->block[id]) += param->delta*adapt*T;
	}
	else
	{
		VectorXd s_(param->A_row);
		VectorXd Ax_(param->A_row);
		int size = data->block[id + 1] - data->block[id];
		int start = data->block[id] - A_row;
		s_ = *s;
		Ax_ = *Ax;

		VectorXd x_(size);
		x_=(*x).segment(start, size);
		VectorXd temp;
		temp = data->c.segment(start, size) + (data->A.transpose()).block(start, 0, size, param->A_row) * (s_ + 2*data->gamma.cwiseProduct(Ax_ - data->b));
		temp = (x_ - data->eta.segment(start, size).cwiseProduct(temp));
		temp = 0.5*temp+0.5*temp.cwiseProduct(temp.cwiseSign());
		(*x).segment(start, size) += param->delta*adapt*(temp-x_);
		std::lock_guard<std::mutex> guard(mu);
		(*Ax) += param->delta*adapt*data->A.block(0, start, A_row, size)*(temp-x_);
	}
}

void PD::parallel_stop(int id)
{
	int worker_num = param->total_num_threads;
	// copy Ax
	int start = (A_row / worker_num)*id;
	int size = id < worker_num - 1 ? A_row / worker_num : A_row - (worker_num - 1)*(A_row / worker_num);
	VectorXd Ax_(size);
	// compute ATs
	start = (A_col / worker_num)*id;
	size = id < worker_num - 1 ? A_col / worker_num : A_col - (worker_num - 1)*(A_col / worker_num);
	//VectorXd temp(size);
	(*ATs).segment(start, size) = ((data->A.block(0, start, A_row, size)).transpose()*(*s))+data->c.segment(start,size);
}

void PD::check_stop()
{
	VectorXd temp(param->A_col);
	temp = (*ATs);
 	temp = temp*0.5 - 0.5*temp.cwiseProduct(temp.cwiseSign());
	crit->dfeas = temp.lpNorm<Eigen::Infinity>();
	crit->pfeas = (*Ax-data->b).lpNorm<Eigen::Infinity>();
	double primal = (*x).dot(data->c);
	double dual = -(*s).dot(data->b);
	crit->relgap = fabs(primal - dual) / (1 + fabs(primal) + fabs(dual));
	if (crit->pfeas<param->feastol*(1+data->b_linf) && crit->dfeas < param->feastol*(1+data->c_linf)&& crit->relgap < param->feastol){
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

void PD::get_minizer()
{
	return;
}
