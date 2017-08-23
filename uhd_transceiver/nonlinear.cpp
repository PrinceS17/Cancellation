#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <complex>
#include "wavetable.hpp"

using namespace std;
using namespace Eigen;

MatrixXf x2A(VectorXf & x, int l, int k, int dim)     //A(0,0) is x(0), 1st y is y(k)
{
	int n = l - 1;				// l: length of [x0,...,xn]
	int st = k;					// start of y, x.size() should >= n + 2k, better =
	//int dim = 4;				// dimension of nonlinearity
	
	MatrixXf A(n + 1, 2*k*dim);
	for(int i = 0; i < n + 1; i ++)
	for(int j = 0; j < 2*k; j ++)
	for(int p = 0; p < dim; p ++)
		A(i, j*dim + p) = pow(x(st - k + j + i),p + 1);		// x^(p + 1)
	return A;
}

MatrixXf x2A(VectorXf & x, int k, int dim)
{
	int lx = x.size();
	return x2A(x, lx - 2*k + 1, k, dim);
}

MatrixXf x2A(VectorXf & x, int k)	        // default dim is 1
{
	return x2A(x, k, 1);
}