#pragma once
#ifndef ALP_INCLUDE_INDEX
#define ALP_INCLUDE_INDEX

#include "algebra.h"
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;
using namespace Eigen;

struct Parameter
{
	// info of iteration
	int max_iter;
	int total_num_threads;
	int block_size_p;
	int block_size_d;
	int block_num;
	int check_step;	// every *this check optiaml

	// status
	bool stop;

	// info of the problem
	int A_col;
	int A_row;

	// parameters
	double eta;
	double delta;
	double alpha;
	double feastol;	// optimal

	// method
	bool scale;	// whether to scale
	int choice;	// 
	// initial value
	Parameter() : eta(1.0), max_iter(5000), block_size_p(1), block_size_d(1), feastol(1e-5),
		scale(false), check_step(100), delta(0.5), alpha(1.0), stop(false) {}
};

struct Data
{
	// primitive data
	MatrixXd A;
	VectorXd b;
	VectorXd c;
	VectorXd ans;
	double b_linf;
	double c_linf;
	// pre-condition
	VectorXd gamma;
	VectorXd eta;
	// block info
	std::vector<int> block;
};

struct Result
{
	int status;	// 0--processing; 1--optimal;
	double optimal;
	double gap;
	double time;
	double time_without;
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
