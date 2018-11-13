#ifndef ALP_INCLUDE_INDEX
#define ALP_INCLUDE_INDEX

#include "algebra.h"
#include <vector>
struct Parameter
{
	// info of iteration
	int max_iter;
	int total_num_threads;
	int block_size;
	int block_num;
	int check_step;	// every *this check optiaml

	// status
	bool stop;

	// info of the problem
	int A_col;
	int A_row;

	// parameters
	double eta;
	double lambda;
	double feastol;	// optimal

	// method
	bool scale;	// whether to scale A
	int choice;	// 
	
	// initial value
	Parameter() : eta(1.0), max_iter(5000), block_size(1), feastol(1e-5),
		scale(false), check_step(100), lambda(1.0), stop(false) {}
};

using namespace Eigen;
struct Data
{
	// primitive data
	MatrixXd A;
	VectorXd b;
	VectorXd c;
	double c_linf;
	double b_linf;

	// intermediate
	VectorXd d;
	MatrixXd M;
	MatrixXd T;
	// block info
	std::vector<int> block;
};

struct Result
{
	int status;	// 0--processing; 1--optimal;
	double optimal;
	double gap;
	double time;
	double check_time;
	int iter;
	VectorXd minimizer;
	Result() : status(0), optimal(0.), gap(0.), time(0.), iter(0) {}
};

struct StoppingCriteria
{
	double pfeas;
	double dfeas;
	double relgap;
	StoppingCriteria() : pfeas(0.), dfeas(0.), relgap(0.) {}
};

#endif
