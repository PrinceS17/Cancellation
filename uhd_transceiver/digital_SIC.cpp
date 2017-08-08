#include <Eigen/Dense>
#include <Eigen/SVD>
#include <vector>
#include "wavetable.hpp"
#include <iostream>
#include <fstream>
#include <complex>
#include <stdint.h>
#include <cmath>
using namespace Eigen;
using namespace std;


// dg_sic begin

// cancellation function: called in main
// 1, matrix A and A_inv: now calculated by function, in fact should store into a global object
// 2, MatrixXf: now use float real type, MatrixXcf probably needed in the future 

MatrixXf x2A(VectorXf& x, int l, int k)		// Note A(0,0) is x(0), so the 1st y is y(k)
{
	int n = l - 1;							// l is the length of [x0, x1, ..., xn]
	int st = k;								// normally x.size() should not be less than n + 2k
	MatrixXf A(n + 1, 2*k);					// should be check a third time!
	for(int i = 0; i < n + 1; i++)
	for(int j = 0; j < 2*k; j++)
	{
		A(i,j) = x(st - k + j + i);
 	}
	return A;
}

MatrixXf x2A(VectorXf& x, int k)
{
	int lx = x.size();	
	MatrixXf A = x2A(x, lx - 2*k + 1, k);			   // every points are used
	return A;
}

Index dg_sync(VectorXf& preamble, VectorXf& rbuff)   // return delay in rbuff
{
	int cor_length = rbuff.size() - preamble.size() + 1;
	VectorXf Cor(cor_length);					   // make the preamble
	for(int i = 0; i < cor_length; i++) 
		Cor(i) = preamble.transpose()*rbuff.segment(i,preamble.size());
	Index idx;								   // simply find the max number
	Cor.maxCoeff(&idx);
	return idx - 1;
}

VectorXf estimate(VectorXf& sbuff, VectorXf& rbuff, int estimator_length)  // 2, estimation: sbuff: -k, ..., n+k-1; rbuff: 0,...,n (now estimate with A calculating repeatly! )
{
	// definition
	int k = estimator_length/2;
	if(sbuff.size() != rbuff.size() + 2*k - 1) 
	{	
		cout<<"\n error: length of pilot and rx signal don't match!"<<endl;
		exit(0);
	}

	// generate A
	MatrixXf A = x2A(sbuff, k);	
	BDCSVD<MatrixXf> svd(A, ComputeThinU|ComputeThinV);		  // use BDC SVD which is for large matrix
	VectorXf h = svd.solve(rbuff);
	
	/* MatrixXf A2 = A.transpose()*A;					   // this way may not work for big estimator length
	if(!A2.determinant())
	{
		cout<<"\n error: A'A is not invertible!"<<endl;
		exit(0);
	}
	h = A2.inverse()*A.transpose()*rbuff;			   // since A's psuedo inverse is (A'A)^-1 * A', it's A_inv*y*/
	
	return h;

}

VectorXf dg_cancel(VectorXf& sbuff, VectorXf& rbuff, VectorXf& h, int estimator_length)  //3, cancellation
{
	// definition
	int k = estimator_length/2;

	// process for buff here

	if(sbuff.size() != rbuff.size() + 2*k - 1) 
	{	
		cout<<"\n error: length of pilot and rx signal don't match!"<<endl;
		exit(0);
	}

	// generate A1
	MatrixXf A1 = x2A(sbuff, k);
	return rbuff - A1*h;

}

VectorXf dg_sic(
	VectorXf &x, 
	VectorXf &y,							       // initial signal got from UHD: here haven't defined complex number
	VectorXf &preamble,					       // should have complete definition later
	int estimator_length,
	int preamble_length,
	int pilot_length,
	int signal_length,						   // the length of TX which made signal_length + delay <= RX's length, may be redefine
	int samp_rate
	)
{
	// print basic information before cancellation
	cout<<endl<<"-----------------start cancellation------------------------------"<<endl;
	cout<<"-- sampling rate: "<<samp_rate<<endl;
	cout<<"-- signal_length: "<<signal_length<<endl;
	cout<<"-- estimator_length: "<<estimator_length<<endl;
	cout<<"-- preamble_length: "<<preamble_length<<endl;
	cout<<"-- pilot_length: "<<pilot_length<<endl<<endl;
	

	int k = estimator_length/2;
	Index delay = dg_sync(preamble, y);				// error here: Index type? or calculate?
	cout<<"-- delay = "<<delay<<endl;

	// define tx&rx_pilot and estimate h
	// error here: Segmentation fault (core dumped)
	
	if(preamble_length + pilot_length + k - 2 >= x.size() | delay + preamble_length + pilot_length - 1 >= y.size())
{
	cout<<"\n error: the last index of requested pilot is beyond given signal!"<<endl;
	exit(0);
}
	VectorXf tx_pilot = x.segment(preamble_length - k, pilot_length + 2*k - 1);
	VectorXf rx_pilot = y.segment(delay + preamble_length, pilot_length);
	VectorXf h = estimate(tx_pilot, rx_pilot, estimator_length);
	//cout<<"h = \n"<<h.transpose()<<endl;


	// obtain new sequence of TX and RX
	int L = signal_length -delay + k -1;		// possible largest length of data for sine wave
	VectorXf tx_data = x.segment(preamble_length + pilot_length - k, L - pilot_length - preamble_length + k);
	VectorXf rx_data = y.segment(delay + preamble_length + pilot_length, L - pilot_length - preamble_length - k + 1);
	VectorXf y_clean;
	y_clean = dg_cancel(tx_data, rx_data, h, estimator_length);

	cout<<"-- x's norm: "<<x.norm()<<endl;
	cout<<"-- y's norm: "<<y.norm()<<endl;	
	cout<<"-- y_clean's norm: "<<y_clean.norm()<<endl;

	// write to file and some other work
	ofstream outfile;
	string name[3] = {"tx_file","rx_file","y_clean_file"};
	float * ptr[3] = {x.data(), y.data(), y_clean.data()};									 // may not work
	for (int i = 0; i < 3; i++)
	{
		outfile.open(name[i].c_str(), ios::out | ios::binary);				 // for Matlab visualization
		if(i < 2)
			outfile.write((const char*)ptr[i], signal_length*sizeof(float));     // try I/Q channel by one first
		else outfile.write((const char*)ptr[i], rx_data.size()*sizeof(float));       // y_clean is shorter than x, y
		outfile.close();
	}

}


