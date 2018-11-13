#ifndef ADRS_INCLUDE_DRS
#define ADRS_INCLUDE_DRS

#include "algebra.h"
#include "index.h"

void parse_input_argv_DRS(Parameter & param, int argc, char* argv[], 
	std::string &A_filename, std::string &B_filename, std::string &c_filename, Data &data);
void pre_computing_DRS(Data&, Parameter&);
void print_parameters(Parameter &param);
class DRS
{
public:
	VectorXd* z;
	VectorXd* s;
	VectorXd* Ax;
	VectorXd* ATs;
	VectorXd* mw;
	Data * data;
	Parameter * param;
	StoppingCriteria * crit;
	Result * res;
	int A_col;
	int A_row;
	int block_num;

	DRS(VectorXd*, VectorXd*, VectorXd*, VectorXd*, VectorXd*, Data*, Parameter*, StoppingCriteria*, Result*);
	~DRS();

	int get_id();
	void update(int);
	void parallel_stop(int id);
	void parallel_getmin(int id);
	void check_stop();
	void get_minizer();
private:
	// temp variable
	VectorXd x;
	VectorXd w;
	VectorXd u;
};


#endif // !ADRS_INCLUDE_DRS
