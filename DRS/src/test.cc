#include <iostream>
#include <thread>
#include <pthread.h>
//#include "pthread.h"
#include <atomic>

#include "util.h"
#include "index.h"
#include "DRS.h"
#include "algebra.h"
#include "worker.h"
using namespace std;
using namespace Eigen;

std::atomic<int> block_id(0);
std::atomic<int> iter(1);
pthread_barrier_t barrier;

int main(int argc, char *argv[])
{
	//cout << "-----program start------";
	// Step 1: Initialization
	Parameter param;
	Data data;
	StoppingCriteria crit;
	Result res;
	std::string A_filename;
	std::string B_filename;
	std::string C_filename;

	parse_input_argv_DRS(param, argc, argv, A_filename, B_filename, C_filename, data);

	LoadMatrix(data.A, A_filename);
	LoadVector(data.b, B_filename);
	LoadVector(data.c, C_filename);
	param.A_col = data.A.cols();
	param.A_row = data.A.rows();

	VectorXd z;
	z = VectorXd::Zero(param.A_col);
	VectorXd s(param.A_row);
	VectorXd Ax(param.A_row);
	VectorXd ATs(param.A_col);
	VectorXd mw(param.A_col);
	pthread_barrier_init(&barrier, NULL, param.total_num_threads);
	cout << "-----------Finish Initialization-----------" << endl;

	// Step 2 Pre-computing
	res.time = get_wall_time();
	pre_computing_DRS(data, param);
	res.time = get_wall_time() - res.time;
	cout << "Pre-computing time is:" << res.time << endl;

	// Step 3 Initialize ADRS
	DRS drs(&z, &s, &Ax, &ATs, &mw, &data, &param, &crit, &res);

	// Step 4 Run
	cout << "---------------------------------" << endl;
	cout << "Calling A-LP: "<< endl;
	cout << "---------------------------------" << endl;

	std::vector<std::thread> threads;
	for (int i = 0; i < param.total_num_threads; i++)
	{
		threads.push_back(std::thread(async_worker, i, drs, std::ref(param), std::ref(res)));
	}
	for (int i = 0; i < param.total_num_threads; i++)
	{
		threads[i].join();
	}
	res.time = get_wall_time()-res.time;

	// Step 5 Print result
	print_parameters(param);
	drs.get_minizer();
	cout << "Computing time is: " << res.time << endl;
	cout << "Computing time without check is:" << res.time - res.check_time << endl;
	cout << "---------------------------------" << endl;
	cout << "iter =: " << iter << endl;
	cout << "status =: " << res.status << endl;
	cout << "optimal =: " << res.optimal << endl;
	cout << "------------Finish---------------" << endl;
	print_vector(z);
	return 0;
}
