#include "algebra.h"
#include <string>
#include <iostream>
#include <fstream>
#include <set>
using namespace std;
using namespace Eigen;

void LoadMatrix(MatrixXd& A, std::string& filename)
{
	ifstream A_file(filename.c_str());
	if (!A_file)
	{
		cout << "cannot open data file" << filename << endl;
	}
	//clear();
	string line;
	getline(A_file, line);
	if (line.substr(0, 15) != "%%MatrixMarket ")
	{
		cout << "header: Not a Matrix Market file" << endl;
	}
	istringstream in(line.substr(15));
	string word;
	set<string> opt;
	while (in)
	{
		in >> word;
		for (unsigned i = 0; i < word.size(); i++)
		{
			word[i] = tolower(word[i]);
		}
		opt.insert(word);
	}
	while (A_file.peek() == '%')
	{
		getline(A_file, line);
	}
	if (opt.find("matrix") == opt.end())
	{
		cout << "header: Not a matrix" << endl;
	}
	if (opt.find("real") == opt.end())
	{
		cout << "header: Not a real matrix" << endl;
	}
	if (opt.find("general") == opt.end())
	{
		cout << "header: Not a general matrix" << endl;
	}
	int n, m;
	A_file >> n >> m;
	A.resize(n, m);
	double val;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A_file >> val;
			A(j, i) = val;
		}
	}
	A_file.close();
}

void LoadVector(VectorXd& A, std::string& filename)
{
	ifstream A_file(filename.c_str());
	if (!A_file)
	{
		cout << "cannot open data file" << filename << endl;
	}
	// ?
	//clear();
	string line;
	getline(A_file, line);
	if (line.substr(0, 15) != "%%MatrixMarket ")
	{
		cout << "header: Not a Matrix Market file" << endl;
	}
	istringstream in(line.substr(15));
	string word;
	set<string> opt;
	while (in)
	{
		in >> word;
		for (unsigned i = 0; i < word.size(); i++)
		{
			word[i] = tolower(word[i]);
		}
		opt.insert(word);
	}
	while (A_file.peek() == '%')
	{
		getline(A_file, line);
	}
	if (opt.find("matrix") == opt.end())
	{
		cout << "header: Not a matrix" << endl;
	}
	if (opt.find("real") == opt.end())
	{
		cout << "header: Not a real matrix" << endl;
	}
	if (opt.find("general") == opt.end())
	{
		cout << "header: Not a general matrix" << endl;
	}
	int n, m;
	A_file >> n >> m;
	A.resize(n);
	double val;
	for (int j = 0; j < n; j++)
	{
		A_file >> val;
		A(j) = val;
	}
	A_file.close();
}

void print_vector(VectorXd & a)
{
	for (int i = 0; i < a.size(); i++)
	{
		cout << a(i) << " ";
	}
	cout << endl;
}

