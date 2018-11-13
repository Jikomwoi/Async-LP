#include <iostream>
#include <thread>
#include <pthread.h>
#include <atomic>
#include <mutex>
#include "index.h"
#include "primdual.h"
#include "algebra.h"
#include "worker.h"
#include "util.h"
using namespace std;
using namespace Eigen;

std::mutex mu;
std::atomic<int> block_id(0);
std::atomic<int> iter(1);
pthread_barrier_t barrier;

int main(int argc, char *argv[])
{
	// Step 1: Initialization
	Parameter param;
	Data data;
	StoppingCriteria crit;
	Result res;
	std::string A_filename;
	std::string B_filename;
	std::string C_filename;
	parse_input_argv(param, argc, argv, A_filename, B_filename, C_filename,  data);
	
	LoadMatrix(data.A, A_filename);
	LoadVector(data.b, B_filename);
	LoadVector(data.c, C_filename);
	param.A_col = data.A.cols();
	param.A_row = data.A.rows();

	VectorXd s;
	s = VectorXd::Zero(param.A_row);
	VectorXd x;
	x = VectorXd::Zero(param.A_col);
	VectorXd Ax;
	Ax = VectorXd::Zero(param.A_row);
	VectorXd ATs;
	ATs = VectorXd::Zero(param.A_col);
	pthread_barrier_init(&barrier, NULL, param.total_num_threads);

	srand((unsigned int)(time(NULL)));
	//cout << "-----------Finish Initialization-----------" << endl;

	// Step 2 Pre-computing
	res.time = get_wall_time();
	pre_computing_PD(data, param);
	res.time = get_wall_time() - res.time;
	cout << "Pre-computing time is:" << res.time << endl;

	// Step 3 Initialize ALP
	PD pd(&s, &x, &Ax, &ATs, &data, &param, &crit, &res);

	// Step 4 Run
	cout << "---------------------------------" << endl;
	cout << "Calling A-LP: " << endl;
	cout << "---------------------------------" << endl;
	std::vector<std::thread> threads;
	for (int i = 0; i < param.total_num_threads; i++)
	{
		threads.push_back(std::thread(async_worker, i, pd, std::ref(param), std::ref(res)));
	}
	for (int i = 0; i < param.total_num_threads; i++)
	{
		threads[i].join();
	}
	res.time = get_wall_time() - res.time;
	
	// Step 5 Print result
	print_parameters(param);
	pd.get_minizer();
	cout << "Computing time is: " << res.time << endl;
	cout << "Computing time without check: " << res.time - res.time_without << endl;
	cout << "---------------------------------" << endl;
	cout << "iter =: " << iter << endl;
	cout << "status =: " << res.status << endl;
	cout << "optimal =: " << res.optimal << endl;
	cout << "---------------------------------" << endl;
	cout << x << endl;
	return 0;
}

