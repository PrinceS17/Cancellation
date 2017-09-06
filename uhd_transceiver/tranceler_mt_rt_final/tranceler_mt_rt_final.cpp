/* ======================tranceler_mt_rt_final.cpp =========================
This is the UHD code of multi-tone sine wave real time transceiver & nonlinear canceler. It is tested with USRP B205mini and based on VectorXf of Eigen. 
Different from UHD examples, it uses VectorXf to generate sine wave of precise frequency rather than use wavetable.hpp.
It depends on fft.hpp, nonlinear_peak.hpp. fft() and sic_db() are in fft.hpp, and nonlinear x2A() and peaks() are in nonlinear_peak.hpp.
*/

#include "nonlinear_peak.hpp"
#include "fft.hpp"
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>
#include <uhd/types/tune_request.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <boost/signals2/mutex.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <vector>
#include <iostream>
#include <fstream>
#include <complex>
#include <csignal>
#include <stdint.h>
#include <cmath>

#define pi 3.14159265358979323846

using namespace Eigen;
using namespace std;
namespace po = boost::program_options;


static bool stop_signal_called = false;		// control signal of all threads
boost::mutex mtx;							// mutex for global variables below
VectorXcf global_tx, global_rx;				// long global buffers for the storage of TX/RX buffer
int global_num[2] = {0, 0};					// # of TX/RX buffer in global_tx, rx
int gsize;

void sig_int_handler (int) {stop_signal_called = true;}

// generate mt_sine to send
VectorXcf mt_sine_generation(VectorXf &w_freq, VectorXf &wave_freq, int wave_num, double wave_freq_1, double wave_space, int spb, double rate, double ampl)
{	

	VectorXf t = VectorXf::LinSpaced(spb, 1, spb)/rate;
	VectorXf temp = VectorXf::Zero(spb);
	if(!w_freq(0))				// if freq-1 is 0, then choose 1st way to generate mt_sine
		 wave_freq = VectorXf::LinSpaced(wave_num, wave_freq_1, wave_freq_1 + (wave_num - 1)*wave_space);
	else 					// else, get the frequency from input
	{	
		wave_num = (w_freq.array() > 0).count();
		wave_freq = w_freq.segment(0, wave_num);	
	}
	for(int i = 0; i < wave_num; i ++)
		temp.array() = temp.array() + ampl / wave_num * (2*pi*wave_freq[i]*t).array().sin() ;
	VectorXcf mt_sine = temp.cast<complex <float> > ();
	return mt_sine;

}

// digital canceler part
int dg_sync(VectorXf &preamble, VectorXf &rbuff)	// use preamble to get delay in rbuff
{
	int cor_length = rbuff.size() - preamble.size() + 1;
	if(cor_length < 0)
	{
		cout<<"-- Error: RX signal is shorter than preamble!"<<endl;
		exit(0);
	}
	VectorXf Cor(cor_length);
	for(int i = 0; i < cor_length; i ++)
		Cor(i) = preamble.transpose() * rbuff.segment(i, preamble.size());
	int peak_num = 7;
	VectorXf idxs = peaks(Cor, peak_num);						// find the biggest peak_num peaks index
	sort(idxs.data(), idxs.data() + peak_num, less<int>());		// find the earliest peak among peaks above
	return idxs(0);
}

VectorXf estimate(VectorXf &sbuff, VectorXf &rbuff, int estimator_length)
{
	int k = estimator_length / 2;
	if(sbuff.size() != rbuff.size() + 2*k - 1)
	{
		cout<<"-- Error: length of pilot and RX signal don't match!"<<endl;
		exit(0);
	}
	MatrixXf A = x2A(sbuff, k);
	BDCSVD <MatrixXf> svd(A, ComputeThinU | ComputeThinV);
	VectorXf h = svd.solve(rbuff);							// LSQ: calculate h (same as using psuedo inverse)
	return h;
}

VectorXf dg_cancel(VectorXf &sbuff, VectorXf &rbuff, VectorXf &h, int estimator_length) 
{
	int k = estimator_length / 2;
	if(sbuff.size() != rbuff.size() + 2*k - 1)
	{
		cout<<"-- Error: length of pilot and RX signal don't match!"<<endl;
		exit(0);
	}
	MatrixXf A1 = x2A(sbuff, k);
	return rbuff - A1*h;
}

VectorXf digital_canceler(
	VectorXf &preamble,	
	VectorXf &wave_freq,
	int samp_rate,
	int spb,
	int start,
	int *len,
	string *file
	)
{	
	// print basic information before cancellation
	int estimator_length = len[0];
	int preamble_length = len[1];
	int pilot_length = len[2];			// notice: pilot length is n + 1
	int signal_length = len[3];	

	cout<<endl<<"------------------- start cancellation -------------------------"<<endl;
	cout<<"-- sampling rate: "<<samp_rate<<endl;
	cout<<"-- samples per buffer: "<<spb<<endl;
	cout<<"-- signal_length: "<<signal_length<<endl;
	cout<<"-- start: "<<start<<endl;
	cout<<"-- estimator_length: "<<estimator_length<<endl;
	cout<<"-- preamble_length: "<<preamble_length<<endl;
	cout<<"-- pilot_length: "<<pilot_length<<endl<<endl;
	
	int delay_const = 1.2e5;		// default constant of TX, RX delay
	int can_num = 0;					// # of cancellation iterations
	int num[3] = {0, 0, 0};			// # of samples written into file
	ofstream outfile[3];
	for(int i = 0; i < 3; i ++)
		outfile[i].open(file[i].c_str(), ios::out | ios::binary);
	VectorXf x(signal_length), y(signal_length);
	float result_av = 0;

	while(! stop_signal_called)
	{
		// read data from global buffers
		mtx.lock();
		int tx_num = global_num[0];
		int rx_num = global_num[1];
		
		if(rx_num <= delay_const / spb + 3)		// avoid reading 0s when there's no signal in global buffers
		{
			mtx.unlock();
			continue;
		}
		int can_pos = can_num * signal_length;
		for(int i = 0; i <signal_length; i ++)	
		{
			x(i) = global_tx[(can_pos + i) % gsize].real();
			y(i) = global_rx[(can_pos + i + delay_const) % gsize].real();
		}
		can_num ++;	
		mtx.unlock();
		
		
		// define TX/RX pilot & estimate h
		int k = estimator_length / 2;
		int delay = dg_sync(preamble, y);
		if(preamble_length + pilot_length + k - 1 > x.size() | delay + preamble_length + pilot_length  > y.size())
		{
			cout<<"-- Error: requested pilot exceeds boundary!"<<endl;
			exit(0);
		}
		if(start + preamble_length < k)
		{
			cout<<"-- Error: no enough samples for channel estimation!"<<endl;
			exit(0);
		}
		VectorXf tx_pilot = x.segment(start + preamble_length - k, pilot_length + 2*k - 1);
		VectorXf rx_pilot = y.segment(start + delay + preamble_length, pilot_length);
		VectorXf h = estimate(tx_pilot, rx_pilot, estimator_length);


		// do the cancellation
		int L = signal_length - delay - start + k -1;		// possible largest length of data for sine wave
		VectorXf tx_data = x.segment(start + preamble_length + pilot_length - k, L - pilot_length - preamble_length + k);
		VectorXf rx_data = y.segment(start + delay + preamble_length + pilot_length, L - pilot_length - preamble_length - k + 1);
		VectorXf y_clean = dg_cancel(tx_data, rx_data, h, estimator_length);


		// calculate the cancellation and std output
		double bw = 2 * samp_rate / (double)signal_length;
		double rg = 5e3;
		int wave_num = wave_freq.size();	
		int prd = 20;
		if(can_num % prd == 0)
		{	
			cout<<"-- TX No. "<<tx_num<<" , RX No. "<<rx_num<<" , Cancel No. "<<can_num * signal_length / spb<<endl;
			
			//MatrixXd ary(wave_num, 2);
			double **ary = new double *[wave_num];
			VectorXf result(wave_num);			
			for(int i = 0; i < wave_num; i ++)
			{
				ary[i] = new double [2];
				result[i] = sic_db(y, y_clean, samp_rate, wave_freq[i], bw, rg, ary[i]);
				cout<<"   "<<wave_freq[i]/1e3<<" kHz: "<<setprecision(3)<< result[i]<<" dB, "<< ary[i][0]<<" -> "<< ary[i][1]<<" dB"<<endl;
			}
			cout<<"   Total: "<<setprecision(3)<<result.mean()<<" dB"<<endl;
			
			// calculate average cancellation result
			if(can_num / prd < 11) result_av += result.mean();
			if(can_num == prd*10) 
			{
				result_av /= 10;
				cout<<"-- result_av == "<<result_av<<endl<<endl;
			}
		}

		
		
		
		// write results to file
		float * ptr[3] = {x.data(), y.data(), y_clean.data()};
		int length[3] = {x.size(), y.size(), y_clean.size()};
		for (int i = 0; i < 3; i++)
		{
			outfile[i].write((const char*)ptr[i], length[i]*sizeof(float));     // try I/Q channel by one first
			num[i] += length[i];

			if(num[i]*sizeof(float) > 100e6*sizeof(char))		// update when size increases to limit
			{					
				outfile[i].close();
				outfile[i].open(file[i].c_str(),ios::binary | ios::out);
				num[i] = 0;
				cout<<"-- New "<<file[i]<<endl;
			}
		}
	}
	for(int i = 0; i < 3; i ++)
		outfile[i].close();
}

void transmitter(
	VectorXcf mt_sine,			// multi-tone sine wave
	uhd::tx_streamer::sptr tx,
	uhd::tx_metadata_t md,
	int spb
	)
{
	int num = 0;
	complex <float> * st = mt_sine.data();
	vector <complex <float> > buff(st, st + mt_sine.size());

	while(! stop_signal_called)
	{
		mtx.lock();
		int tx_pos = spb*global_num[0];
		for(int i = 0; i < spb; i ++)
			global_tx[(tx_pos + i) % gsize] = buff[i];
		global_num[0] ++;
		mtx.unlock();

		num += tx->send(&buff.front(), buff.size(), md);		// send TX data from TX buffer
		md.start_of_burst = false;
		md.has_time_spec = false;
	}

	// end transmitting
	md.end_of_burst = true;
	tx->send("", 0, md);
	mtx.lock();
	cout<<"-- TX done!"<<endl;
	mtx.unlock();
}

void receiver(
	uhd::rx_streamer::sptr rx,
	int spb,
	double num_requested_samples,
	double time_requested = 0,
	bool enable_size_map = false,
	bool continue_on_bad_packet = false
	)
{
	uhd::rx_metadata_t md;
	vector<complex <float> > rbuff(spb);
	int num_total_samps = 0;

	// set up RX streaming
	uhd::stream_cmd_t stream_cmd((num_requested_samples == 0)?
		uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS:
	uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE
		);
	stream_cmd.num_samps = size_t(num_requested_samples);
	stream_cmd.stream_now = true;
	stream_cmd.time_spec = uhd::time_spec_t();
	rx->issue_stream_cmd(stream_cmd);                // tells the usrp to send samples to host

	typedef std::map<size_t,size_t> SizeMap;
	SizeMap mapSizes;
	bool overflow_message = true;

	while(! stop_signal_called && (num_requested_samples!= num_total_samps || !num_requested_samples))
	{
		int num_rx_samps = rx->recv(&rbuff.front(), rbuff.size(), md, 0.1f);
		mtx.lock();
		int rx_pos = spb*global_num[1];
		for(int i = 0; i < spb; i ++)
			global_rx[(rx_pos + i) % gsize] = rbuff[i];
		global_num[1] ++;
		mtx.unlock();

		if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) 
		{
			std::cout << boost::format("Timeout while streaming") << std::endl;
			break;
		}
		if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW)
		{
			if (overflow_message) 
			{
				overflow_message = false;
				std::cerr << boost::format(
					"Got an overflow indication."
					) ;
			}
			continue;
		}
		if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE)
		{
			std::string error = str(boost::format("Receiver error: %s") % md.strerror());
			if (continue_on_bad_packet)
			{
				std::cerr << error << std::endl;
				continue;
			}
			else
				throw std::runtime_error(error);
		}		// errors
		if (enable_size_map) 
		{
			SizeMap::iterator it = mapSizes.find(num_rx_samps);
			if (it == mapSizes.end())
				mapSizes[num_rx_samps] = 0;
			mapSizes[num_rx_samps] += 1;
		} 		
		num_total_samps += num_rx_samps;
	}
	// end receiving
	uhd::stream_cmd_t rx_cmd(uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS);
	rx->issue_stream_cmd(rx_cmd);
	mtx.lock();
	cout<<"-- RX done!"<<endl;
	mtx.unlock();
}


int UHD_SAFE_MAIN(int argc, char *argv[])
{
	string cpu, wire, subdev, a1, a2, ref, pps, tx_args, rx_args;
	cpu = "fc32";
	wire = "sc16";
	subdev = "A:A";
	a1 = "TX/RX";
	a2 = "RX2";
	ref = "internal";
	pps = "internal";
	string file[3] = {"tx_file", "rx_file", "y_clean_file"};

	double freq, gain, total_num_samps, total_time, bw;		// default setting 
	freq = 915e6;
	gain = 65;
	bw = 1e6;
	total_num_samps = 0;
	total_time = 20;

	double rate, tx_gain, rx_gain;		// parameters set by po
	float ampl;

	double wave_freq_1, wave_space;								// about multi-sine generation
	int wave_num;										// 4 tones by default
	
	VectorXf w_freq(8);									// for all frequency: should be different from wave_freq

	uhd::set_thread_priority();

	po::options_description desc("Allowed options");
	desc.add_options()
		("tx_args",po::value<string>(&tx_args)->default_value(""),"uhd device address args")
		("rx_args",po::value<string>(&rx_args)->default_value(""),"uhd device address args")
		("rate", po::value<double>(&rate)->default_value(2e6), "rate of transmit and receive samples")
        	("ampl", po::value<float>(&ampl)->default_value(float(0.3)), "amplitude of the waveform [0 to 0.7]")
		("tx-gain", po::value<double>(&tx_gain)->default_value(gain), "gain for the transmit RF chain")
		("rx-gain", po::value<double>(&rx_gain)->default_value(gain), "gain for the receive RF chain")
		("wave-num", po::value<int>(&wave_num)->default_value(4), "number of sine wave tones")
		("wave-freq-1", po::value<double>(&wave_freq_1)->default_value(100e3), "1st waveform frequency in Hz")
		("wave-space", po::value<double>(&wave_space)->default_value(100e3), "spacing between adjacent tones")
		("freq-1", po::value<float>(&w_freq(0))->default_value(0), "1st tone of sine wave")
		("freq-2", po::value<float>(&w_freq(1))->default_value(0), "2nd tone of sine wave")
		("freq-3", po::value<float>(&w_freq(2))->default_value(0), "3rd tone of sine wave")
		("freq-4", po::value<float>(&w_freq(3))->default_value(0), "4th tone of sine wave")
		("freq-5", po::value<float>(&w_freq(4))->default_value(0), "5th tone of sine wave")
		("freq-6", po::value<float>(&w_freq(5))->default_value(0), "6th tone of sine wave")
		("freq-7", po::value<float>(&w_freq(6))->default_value(0), "7th tone of sine wave")
		("freq-8", po::value<float>(&w_freq(7))->default_value(0), "8th tone of sine wave")
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);


	// set the usrp
	uhd::usrp::multi_usrp::sptr tx_usrp = uhd::usrp::multi_usrp::make(tx_args);
	uhd::usrp::multi_usrp::sptr rx_usrp = uhd::usrp::multi_usrp::make(rx_args);

	tx_usrp->set_clock_source(ref);
	tx_usrp->set_tx_rate(rate);
	rx_usrp->set_rx_rate(rate);
	cout<<boost::format("Actual tx rate: %f Msps ...")%(tx_usrp->get_tx_rate()/1e6)<<endl;
	cout<<boost::format("Actual rx rate: %f Msps ...")%(rx_usrp->get_rx_rate()/1e6)<<endl;

	tx_usrp->set_tx_antenna(a1);
	rx_usrp->set_rx_antenna(a2);
	
	tx_usrp->set_tx_bandwidth(bw);
	rx_usrp->set_rx_bandwidth(bw);

	tx_usrp->set_tx_freq(freq);
	rx_usrp->set_rx_freq(freq);

	tx_usrp->set_tx_gain(tx_gain);
	rx_usrp->set_rx_gain(rx_gain);
	
	tx_usrp->set_tx_subdev_spec(subdev);
	rx_usrp->set_rx_subdev_spec(subdev);

	std::cout<<boost::format("Using TX Device: %s") % tx_usrp->get_pp_string()<<endl;
	std::cout<<boost::format("Using RX Device: %s") % rx_usrp->get_pp_string()<<endl;
    

	// set transmit & receive streamer
	uhd::stream_args_t stream_args(cpu, wire);
	uhd::tx_streamer::sptr tx = tx_usrp->get_tx_stream(stream_args);
	uhd::rx_streamer::sptr rx = rx_usrp->get_rx_stream(stream_args);
	int spb = tx->get_max_num_samps() * 10;
	tx_usrp->set_time_now(0.0);
	

	// calculate multi-tone sine wave
	VectorXf wave_freq;
	VectorXcf mt_sine = mt_sine_generation(w_freq, wave_freq, wave_num, wave_freq_1, wave_space, spb, rate, ampl);


	// set global buffers
	int num_buff = 1000;
	global_rx = VectorXcf::Zero(num_buff * spb);
	global_tx = VectorXcf::Zero(num_buff * spb);
	gsize = num_buff * spb;


	// check rdf and L0 lock detect is here
    std::vector<std::string> tx_sensor_names, rx_sensor_names;
    tx_sensor_names = tx_usrp->get_tx_sensor_names(0);
    if (std::find(tx_sensor_names.begin(), tx_sensor_names.end(), "lo_locked") != tx_sensor_names.end()) {
        uhd::sensor_value_t lo_locked = tx_usrp->get_tx_sensor("lo_locked",0);
        std::cout << boost::format("Checking TX: %s ...") % lo_locked.to_pp_string() << std::endl;
        UHD_ASSERT_THROW(lo_locked.to_bool());
    }
    rx_sensor_names = rx_usrp->get_rx_sensor_names(0);
    if (std::find(rx_sensor_names.begin(), rx_sensor_names.end(), "lo_locked") != rx_sensor_names.end()) {
        uhd::sensor_value_t lo_locked = rx_usrp->get_rx_sensor("lo_locked",0);
        std::cout << boost::format("Checking RX: %s ...") % lo_locked.to_pp_string() << std::endl;
        UHD_ASSERT_THROW(lo_locked.to_bool());
    }
    tx_sensor_names = tx_usrp->get_mboard_sensor_names(0);

	if (total_num_samps == 0){
		std::signal(SIGINT, &sig_int_handler);
		std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
	}

	
    //reset usrp time to prepare for transmit/receive
    std::cout << boost::format("Setting device timestamp to 0...") << std::endl;
    tx_usrp->set_time_now(uhd::time_spec_t(0.0));


	// set up metadata for transmitter
	uhd::tx_metadata_t md;
	md.start_of_burst = true;				// set start of burst to true for 1st packet in the chain
	md.end_of_burst = false;                // set end of burst to true for the last packet in the chain
	md.has_time_spec = true;                // set true to send at the time specified by time spec; set false to send immediately
	md.time_spec = uhd::time_spec_t(0.1);
	

	// set up signal structure and generate preamble from multi-tone frequency
	int start = 30;
	int preamble_length = rate/wave_freq_1;
	int estimator_length = 42;
	int pilot_length = 600;
	int signal_length = spb;
	int len[4] = {estimator_length, preamble_length, pilot_length, signal_length};
	VectorXf preamble = mt_sine.segment(0, preamble_length).real();

	
	// set up threads for transmitter and digital canceler
	boost::thread_group threads;
	threads.create_thread(boost::bind(&transmitter, mt_sine, tx, md, spb));
	threads.create_thread(boost::bind(&receiver, rx, spb, total_num_samps, 0, false, false));
	digital_canceler(preamble, wave_freq, rate, spb, start, len, file);


	//finished
	stop_signal_called = true;	
	cout<<"-- Done!"<<endl<<endl;
	return EXIT_SUCCESS;

}
