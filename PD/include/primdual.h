#ifndef ADRS_INCLUDE_PD
#define ADRS_INCLUDE_PD

#include "algebra.h"
#include "index.h"

void parse_input_argv(Parameter & param, int argc, char* argv[],
	std::string &A_filename, std::string &B_filename, std::string &c_filename, Data &data);
void pre_computing_PD(Data&, Parameter&);
void print_parameters(Parameter &param);

class PD
{
public:
	VectorXd * s;
	VectorXd* x;
	VectorXd* Ax;
	VectorXd* ATs;
	Data * data;
	Parameter * param;
	StoppingCriteria * crit;
	Result * res;
	int A_col;
	int A_row;
	int block_num;
	PD(VectorXd*, VectorXd*, VectorXd*, VectorXd*, Data*, Parameter*, StoppingCriteria*, Result*);

	int get_id();
	void update(int, int);
	void parallel_stop(int id);
	void check_stop();
	void get_minizer();

private:
	// temp variable
};


#endif ADRS_INCLUDE_PD
