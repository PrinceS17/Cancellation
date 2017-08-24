#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
using namespace Eigen;
using namespace std;

double sic_db(VectorXf &y, VectorXf &y_clean, double rate, double fc, double bw) 	// fc: center frequency; bw: bandwidth; not tested
{
	int L = max(y.size(), y_clean.size());
	int N_fft = pow(2, (int)log2(L) + 1);
	MatrixXf signal(N_fft, 2);
	// complete the  col vector by 0s for FFT, very stupid code
	signal.col(1).segment(0, y.size()) = y ;
	signal.col(1).segment(y.size(), N_fft - y.size()) = VectorXf::zero(N_fft - y.size()) ;
	signal.col(2).segment(0, y_clean.size()) = y_clean ;
	signal.col(2).segment(y_clean.size(), N_fft - y_clean.size()) = VectorXf::zero(N_fft - y.size()) ;
	
	float P_db[2];

	for(int i = 0; i < 2; i ++)
	{
		VectorXf x = signal.col(i);
		VectorXcf fx(N_fft);
		fft(x.data(), NULL, N_fft, fx.real().data(), fx.imag().data());
		
		VectorXf Px(N_fft), temp(N_fft);
		Px.array() = 10*log10(fx.array().abs2() /100 * 1000);		// calculate power spectrum in dB
		temp.segment(0, N_fft/2) = Px(N_fft/2, N_fft/2);		// fftshift
		temp.segment(N_fft/2, N_fft/2) = Px(0, N_fft/2);
		Px = temp;

		int fl_id = ( (fc - bw/2) /rate + 0.5 )*N_fft;
		int fr_id = ( (fc + bw/2) /rate + 0.5 )*N_fft;
		P_db[i] = Px.segment(fl_id, fr_id - fl_id + 1).mean(); 

	} 
		
	return P_db[0] - P_db[1];


}
