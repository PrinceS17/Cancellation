#include <iostream>
#include <string>
#include <complex>
#include <fstream>
#include <vector>
#include <csignal>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>
#include <uhd/types/tune_request.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <boost/algorithm/string.hpp>
#include <wavetable.hpp>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
namespace po = boost::program_options;

// cancellation without receiving all the packets
// now the sbuff, rbuff are the matching packets after capturing the bound of symbol
// still need read data from file (matrix A and A_inv)

// 1, capture the bound of symbol (/packet?) i.e. relation between TX & RX
// 2, from pilot: calculate pilot information A, A_inv and h, 1) h = A_inv*y
// 3, from signal: use 1,2 to do the cancellation : 2) y' = A'h 

 
MatrixXf x2A(VectorXf x, int l, int k)
{
	int n = l - 1;
	int st = k;								// normally x.size() should not be less than n + 2k
	MatrixXf A(n + 1, 2*k);					// should be check a third time!
	for(int i = 0; i < n + 1; i++)
	for(int j = 0; j < 2*k; j++)
	{
		A(i,j) = x[st - k + j + i];
 	}
	return A;
}

MatrixXf x2A(VectorXf x, int k)
{
	int lx = x.size();	
	return x2A(x,lx - 2*k,k);			   // every points are used
}

// 1, sync ?


//  send following part as an argument
// 	samp_type ** A = new samp_type *[l];
//  for(int i = 0; i < l; i++)
//  A[i] = new samp_type[2*k];

VectorXf estimate(VectorXf sbuff, VectorXf rbuff, int estimator_length)  // 2, estimation
{
	// definition
	int k = estimator_length/2;

	// process for sbuff here

	if(sbuff.size() != rbuff.size()) 
	{	
		cout<<"error: length of pilot and rx signal don't equal!"<<endl;
		exit(0);
	}

	// generate A
	MatrixXf A = x2A(sbuff, k);
	VectorXf h;
	h = (A.transpose()*A).inverse()*A.transpose()*rbuff;	    // since A's psuedo inverse is (A'A)^-1 * A', it's A_inv*y
	return h;

}


VectorXf dg_cancel(VectorXf sbuff, VectorXf rbuff, VectorXf h, int estimator_length)  //3, cancellation
{
	// definition
	int k = estimator_length/2;

	// process for buff here

	if(sbuff.size() != rbuff.size()) 
	{	
		cout<<"error: length of pilot and rx signal don't equal!"<<endl;
		exit(0);
	}

	// generate A1
	MatrixXf A1 = x2A(sbuff, k);
	VectorXf y_clean;
	y_clean = A1*h;
	return y_clean;

}

void main()
{
	// many preparations
	int estimator_length;
	int pilot_length;
	int samp_rate;

	
	VectorXf tx_pilot, rx_pilot;     // should have complete definition later

	// estimate
	VectorXf h;
	h = estimate(tx_pilot, rx_pilot, estimator_length);

	// obtain new sequence of TX and RX
	VectorXf tx_data, rx_data;
	VectorXf y_clean;
	y_clean = dg_cancel(tx_data, rx_data, h, estimator_length);

	// write to file and some other work

}