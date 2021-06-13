#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H
#include "sparse_matrix.h"
#include <vector>
namespace stable_fluids{
//CG法Denseバージョン
void conjugate_gradient(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x, int n);
//CG法sparseバージョン
void conjugate_gradient(const sparse_matrix &A, const std::vector<double> &b, std::vector<double> &x, int n, int max_itr, double eps);
//gauss_seidel法sparseバージョン
void gauss_seidel(const sparse_matrix_with_diagonal_element &A, const std::vector<double> &b, std::vector<double> &x, int n, int max_itr);
}//namespace stable_fluids
#endif
