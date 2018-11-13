#pragma once

#ifndef ALP_INCLUDE_ALGEBRA
#define ALP_INCLUDE_ALGEBRA
#define EIGEN_USE_BLAS
//#define EIGEN_USE_MKL_ALL
#include <Eigen/Dense>

using namespace Eigen;

void LoadMatrix(MatrixXd &, std::string &);
void LoadVector(VectorXd &, std::string &);

void print_vector(VectorXd &);

#endif // !ALP_INCLUDE_ALGEBRA
