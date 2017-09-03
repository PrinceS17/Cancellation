#include <Eigen/Dense>
#include <iostream>
#include <cmath>

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
	return x2A(x, k, 4);
}

class greater1
{
	VectorXf v;
	public:
		greater1(VectorXf v0): v(v0) {};
	bool operator() (int a, int b) const
	{
		return v(a) > v(b);
	}
};


VectorXf peaks(VectorXf& v, int num)	// return the biggest peaks idx; num: number of peaks
{
	
	VectorXf p(v.size());
	int total_num = 0;
	for(int i = 1; i < v.size() - 1; i++)
	{
		if(v[i] > v[i - 1] & v[i] < v[i + 1])
		{
			p[total_num ++] = v[i];
		}
	}	
	VectorXf peaks = p.segment(0,total_num);
	VectorXf idx = VectorXf::LinSpaced(total_num,0,total_num - 1);
	greater1 grt(peaks);
	sort(idx.data(),idx.data() + total_num,grt);

	if(num > total_num)
	{
		cout<<"-- error: requested peaks is more than actual peaks!"<<endl;
		cout<<"-- num = "<<num<<" , total_num = "<<total_num<<endl;
	}
	return idx.segment(0,num);

}
