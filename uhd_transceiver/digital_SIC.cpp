#include <iostream>
#include <string>
#include <complex>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

// cancellation without receiving all the packets
// now the sbuff, rbuff are the matching packets after capturing the bound of symbol
// still need read data from file (matrix A and A_inv)

// 1, capture the bound of symbol (/packet?) i.e. relation between TX & RX						 need to be done
// 2, from pilot: calculate pilot information A, A_inv and h, 1) h = A_inv*y                     done
// 3, from signal: use 1,2 to do the cancellation : 2) y' = A'h                                  done

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

int dg_sync(VectorXf& preamble, VectorXf& rbuff)   // return delay in rbuff
{
	int cor_length = rbuff.size() - preamble.size() + 1;
	VectorXf Cor(cor_length);					   // make the preamble
	for(int i = 0; i < cor_length; i++)
		Cor(i) = preamble.transpose()*rbuff;
	Index *index;								   // simply find the max number
	Cor.maxCoeff(index);
	return *index - 1;
}

VectorXf estimate(VectorXf& sbuff, VectorXf& rbuff, int estimator_length)  // 2, estimation: sbuff: -k, ..., n+k-1; rbuff: 0,...,n
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
	VectorXf h;
	MatrixXf A2 = A.transpose()*A;					   // this way may not be efficient but work, need improving
	if(!A2.determinant())
	{
		cout<<"\n error: A'A is not invertiable!"<<endl;
		exit(0);
	}
	h = A2.inverse()*A.transpose()*rbuff;			   // since A's psuedo inverse is (A'A)^-1 * A', it's A_inv*y
	return h;

}

VectorXf dg_cancel(VectorXf& sbuff, VectorXf& rbuff, VectorXf& h, int estimator_length)  //3, cancellation
{
	// definition
	int k = estimator_length/2;

	// process for buff here

	if(sbuff.size() != rbuff.size() + 2*k - 1) 
	{	
		cout<<"error: length of pilot and rx signal don't match!"<<endl;
		exit(0);
	}

	// generate A1
	MatrixXf A1 = x2A(sbuff, k);
	return A1*h;

}

void main()
{
	// many preparations
	int estimator_length;
	int preamble_length;
	int pilot_length;
	int signal_length;						   // the length of TX which made signal_length + delay <= RX's length, may be redefine
	int samp_rate;
	int k = estimator_length/2

	VectorXf x, y;							   // initial signal got from UHD: here haven't defined complex number
	VectorXf preamble;					       // should have complete definition later

	int delay = dg_sync(preamble, y);

	// define tx&rx_pilot and estimate h
	VectorXf tx_pilot = VectorXf::block(preamble_length - k - 1, 0, pilot_length + 2*k - 1, 0);
	VectorXf rx_pilot = VectorXf::block(delay + preamble_length, 0, pilot_length, 0);
	VectorXf h = estimate(tx_pilot, rx_pilot, estimator_length);

	// obtain new sequence of TX and RX
	VectorXf tx_data = VectorXf::block(preamble_length + pilot_length - k, 0, signal_length - pilot_length - preamble_length + k, 0);
	VectorXf rx_data = VectorXf::block(delay + preamble_length + pilot_length, 0, signal_length - pilot_length - preamble_length - k + 1, 0);
	VectorXf y_clean;
	y_clean = dg_cancel(tx_data, rx_data, h, estimator_length);

	// write to file and some other work
	ofstream outfile;
	string name[3] = {"tx_file","rx_file","y_clean_file"};
	VectorXf * ptr[3] = {&x, &y, &y_clean};									 // may not work
	for (int i = 0; i < 3; i++)
	{
		outfile.open(name[i].c_str(), ios::out | ios::binary);				 // for Matlab visualization
		outfile.write((const char*)ptr[i], signal_length*sizeof(float));     // try I/Q channel by one first
		outfile.close();
	}

}



//void main()
//{
//	// only for test
//	VectorXf x(10),y(7);
//	x<<1, 2, 3, 4, 5, 6, 7, 8, 9, 10;
//	y<<-3, 8, -1, 0, 1, 2, 3;
//	cout<<"x: \n"<<x.transpose()<<endl;
//	cout<<"y: \n"<<y.transpose()<<endl;
//	int k = 2;
//	VectorXf h = estimate(x, y, 2*k);
//	cout<<"h : \n"<<h<<endl;
//	//MatrixXf m1 = x2A(x, x.size()-2*k, k);
//	//cout<<"x to A: \n"<<m1<<endl;
//
//
//
//}