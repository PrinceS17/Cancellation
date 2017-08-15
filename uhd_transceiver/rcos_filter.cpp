#include <vector>
#include <Eigen/Dense>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;
using namespace Eigen;
#define pi (float)3.14159265358979323846

VectorXf rcos_filter(float beta, int sps, int span)			// tested: right
{
	int N = sps*span;		// total number of samples
	if(N%2) 
	{
		cout<<"-- rcf: Total number of samples is odd!"<<endl;
		return VectorXf::Zero(N + 1);
	}
	int delay = N/2;
	VectorXf t(N + 1);
	t = VectorXf::LinSpaced(N + 1,-delay,delay)/(float)sps;
	//VectorXcf tx = t.array().sin()*(pi*beta*t.array()/(1 - (2*beta*t.array()).square())).cos()/(pi*sps*t.array());    // normal range
	
	VectorXf tx(N + 1);
	double eps = 1.49e-8;
	for (int i = 0; i < N + 1; i ++)					// use loop directly
	{
		if(abs(1 - std::pow(2*beta*t(i),2)) > eps && abs(t(i)) > eps)
		{
			tx(i) = sin(pi*t(i))*cos(pi*beta*t(i));
			tx(i)/= (1 - std::pow(2*beta*t(i),2))*(pi*t(i)*(float)sps);
		}
		else if(abs(t(i)) > eps)
			tx(i) = sin(pi/(2*beta))*beta/(2*(float)sps);		// avoid singular value
		else
			tx(i) = (float)1/(float)sps;
	}
	tx = tx/sqrt(tx.squaredNorm());
	return tx;
}


VectorXcf wave_generation(VectorXcf x, float beta, int sps, int span)		// generate waveform for complex 
{

	VectorXcf x_ups = VectorXcf::Zero(x.size()*sps);
	for(int i = 0; i < x.size(); i ++)
		x_ups(i*sps + 1) = x(i);								// upsampling
	VectorXcf rcf = rcos_filter(beta,sps,span);					// generate raised cosine filter
	VectorXcf tx_wave(x_ups.size());
	for(int i = 0; i < x_ups.size(); i ++)
	{
		if(i < rcf.size() - 1)
			tx_wave(i) = rcf.segment(rcf.size() - i - 1,i + 1).transpose()*x_ups.segment(0,i + 1);
		else
			tx_wave(i) = rcf.transpose()*x_ups.segment(i - rcf.size() + 1,rcf.size());			// convolution: since rcf is symmetric, it is equal to correlation
	}
	return tx_wave;

}

void main()
{
	int signal_length = 50;
	int preamble_length = 10;
	int pilot_length = 10;
	int data_length = signal_length - pilot_length - preamble_length;

	float beta = 0.5;
	int sps = 8;
	int span = 4;

	VectorXcf preamble(preamble_length);						// assign 0,1,0,1... to preamble
	preamble.real() = VectorXf::Zero(preamble_length);
	preamble.imag() = VectorXf::Zero(preamble_length);
	
	for(int i = 0; i < preamble_length/2; i ++)
		preamble[2*i] = complex<float>(1,1);

	VectorXcf sig(pilot_length + data_length);
	for(int i = 0; i< pilot_length + data_length; i ++)
		sig(i) = (double)rand()/RAND_MAX > 0.5? 1.0:0.0;
	sig.imag() = sig.real();
	VectorXcf pilot = sig.segment(0,pilot_length);
	VectorXcf x_bit = sig.segment(pilot_length,data_length);

	VectorXcf x(signal_length);
	x.segment(0, preamble_length) = preamble;
	x.segment(preamble_length, pilot_length) = pilot;
	x.segment(preamble_length + pilot_length, data_length) = x_bit;

	VectorXcf tx = wave_generation(x,beta,sps,span);
	cout<<"tx wave preamble: \n"<<tx.segment(0,preamble_length*sps)<<endl;
}