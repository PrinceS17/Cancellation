#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>
#include <uhd/types/tune_request.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <vector>
#include "wavetable.hpp"
//#include <WinBase.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <csignal>
#include <stdint.h>
#include <cmath>
using namespace Eigen;
using namespace std;
namespace po = boost::program_options;

static bool stop_signal_called = false;
vector<complex<float> > rbuff;
vector<complex<float> > buff;

void sig_int_handler(int){stop_signal_called = true;} 

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



// tranceiver begin

void transmitter(
    //std::vector <complex<float> > &buff,	// no need, transceiver can also generate buff inside the function
    wave_table_class wave_table,
    uhd::tx_streamer::sptr tx0,
    uhd::tx_metadata_t md,
const string file,
    size_t step,
    size_t index,
    int num_channels
)
{

	int num = 0;
	
	// write to file
	ofstream tx_out;
	tx_out.open(file.c_str(),ios::out | ios::binary);
	
	while(not stop_signal_called)	
	{
		for(size_t n = 0; n < buff.size(); n ++)
		{	buff[n] = wave_table(index += step);		  // no problem 
			//buff[n] = n%3;
		}
		num += tx0->send(&buff.front(),buff.size(),md);   // send tx data from buff
		md.start_of_burst = false;
		md.has_time_spec = false;
		
		// 1, output to file
		if (tx_out.is_open())
		{
            		tx_out.write((const char*)&buff.front(), buff.size()*sizeof(float)*2);
			if(num*sizeof(float)*2 > 30e6*sizeof(char)) 
				{						
					tx_out.close();
					tx_out.open(file.c_str(),ios::binary | ios::out);
					num = 0;
					cout<<"New tx_out file!"<<endl;
				}
		}

		// 2, store in global variable



	}
	tx_out.close();	

	// end sending
	
	md.end_of_burst = true;
	tx0->send("",0,md);
	cout<<endl<<"TX done!"<<endl<<endl;
}

template<typename samp_type> void recv_to_file(
	uhd::usrp::multi_usrp::sptr usrp,
    //    uhd::rx_streamer::sptr rx0,
	const string cpu,
	const string wire,
	const string file,
	size_t samps_per_buff,
	double num_requested_samples,
	double time_requested =0,
	bool bw_summary = false,
	bool stats = false,
	bool null = false,
	bool enable_size_map = false,
	bool continue_on_bad_packet = false
	)
{// use cpu and wire type to define the stream
	uhd::stream_args_t stream_args(cpu,wire);
	uhd::rx_streamer::sptr rx0 = usrp->get_rx_stream(stream_args);
	uhd::rx_metadata_t md;
	// vector<samp_type> rbuff(samps_per_buff);    
	std::ofstream outfile;
	outfile.open(file.c_str(),ios::binary | ios::out);
	//string output;
	unsigned long long num_total_samps = 0;

	//set up rx streaming
	uhd::stream_cmd_t stream_cmd((num_requested_samples == 0)?
		uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS:
	uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE
		);
	stream_cmd.num_samps = size_t(num_requested_samples);
	stream_cmd.stream_now = true;
	stream_cmd.time_spec = uhd::time_spec_t();
	rx0->issue_stream_cmd(stream_cmd);                // tells the usrp to send samples to host

	boost::system_time start = boost::get_system_time();
	unsigned long long ticks_requested = (long)(time_requested*(double)boost::posix_time::time_duration::ticks_per_second());
	// it's the whole number of ticks, why?
	boost::system_time last_update = start;
	unsigned long long last_update_samps = 0;
	typedef std::map<size_t,size_t> SizeMap;
	SizeMap mapSizes;
	bool overflow_message = true;


	while(! stop_signal_called && (num_requested_samples!= num_total_samps || !num_requested_samples)){ //! 
		boost::system_time now = boost::get_system_time();
		size_t num_rx_samps = rx0->recv(&rbuff.front(),rbuff.size(),md,0.1f);
		//rx0->recv()

			if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
				std::cout << boost::format("Timeout while streaming") << std::endl;
				break;
			}
			if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW){
				if (overflow_message) {
					overflow_message = false;
					std::cerr << boost::format(
						"Got an overflow indication."
						) ;
				}
				continue;
			}
			if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
				std::string error = str(boost::format("Receiver error: %s") % md.strerror());
				if (continue_on_bad_packet){
					std::cerr << error << std::endl;
					continue;
				}
				else
					throw std::runtime_error(error);
			}  // 3 kinds of errors

			if (enable_size_map) {
				SizeMap::iterator it = mapSizes.find(num_rx_samps);
				if (it == mapSizes.end())
					mapSizes[num_rx_samps] = 0;
				mapSizes[num_rx_samps] += 1;    //confused: what's the original value of mapSizes?
			} 
			num_total_samps += num_rx_samps;
			if (outfile.is_open())
			{
            			outfile.write((const char*)&rbuff.front(), num_rx_samps*sizeof(samp_type));
				if(num_total_samps*sizeof(samp_type) > 30e6*sizeof(char)) 
					{
						outfile.close();
						outfile.open(file.c_str(),ios::binary | ios::out);
						num_total_samps = 0;
cout<<"New output file!"<<endl;
					}
			}
					


	}

// end receiveing
if(outfile.is_open())
	outfile.close();
uhd::stream_cmd_t rx_cmd(uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS);
rx0->issue_stream_cmd(rx_cmd);

}

int UHD_SAFE_MAIN(int argc,char *argv[]){
	// variable definition
	string cpu = "fc32";
	string wire = "sc16";
	string subdev = "A:A";
	string a1 = "TX/RX";
	string a2 = "RX2";
	string ref = "internal";
	string pps = "internal";    // what's pps?
	string tx_args,rx_args,channel_list;
	string type = "float";
	string file = "sine_wave";

	double rate,rx_rate,tx_rate,freq,wave_freq,gain,bw,ampl;
	rate = 2e6;
		file = file +"_" + boost::lexical_cast<string>(rate/1e6) + "M";
	wave_freq = 100e3;
	freq = 915e6;
	gain = 25;				// loop cable: 25; w/o cable: >35
	bw = 1e6;
	rx_rate = rate;
	tx_rate = rate;
	ampl = 0.7;
	double total_num_samps = 0;    // control the total number of samples
	double total_time = 20;

	uhd::set_thread_priority();
	po::options_description desc("Allowed options");
	desc.add_options()
		("help","help message")
		("tx_args",po::value<string>(&tx_args)->default_value(""),"uhd device address args")
		("rx_args",po::value<string>(&rx_args)->default_value(""),"uhd device address args")
		("channels",po::value<string>(&channel_list)->default_value("0"),"which channels to use");


	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

	uhd::usrp::multi_usrp::sptr tx_usrp = uhd::usrp::multi_usrp::make(tx_args);
	uhd::usrp::multi_usrp::sptr rx_usrp = uhd::usrp::multi_usrp::make(rx_args);

	vector<string> channel_strings;
	vector<size_t> channel_nums;
	boost::split(channel_strings, channel_list, boost::is_any_of("\"',")); 	// use " ' or ,to split the string 
	for(size_t ch = 0; ch < channel_strings.size(); ch++){
        size_t chan = boost::lexical_cast<int>(channel_strings[ch]);
        if(chan >= tx_usrp->get_tx_num_channels())
            throw std::runtime_error("Invalid channel(s) specified.");
        else
            channel_nums.push_back(boost::lexical_cast<int>(channel_strings[ch]));  //convert strings into int
    }

	tx_usrp->set_clock_source(ref);
	tx_usrp->set_tx_rate(tx_rate);
	rx_usrp->set_rx_rate(rx_rate);
	cout<<boost::format("Actual tx rate: %f Msps ...")%(tx_usrp->get_tx_rate()/1e6)<<endl;
	cout<<boost::format("Actual rx rate: %f Msps ...")%(rx_usrp->get_rx_rate()/1e6)<<endl;

	tx_usrp->set_tx_antenna(a1);
	rx_usrp->set_rx_antenna(a2);
	
	tx_usrp->set_tx_bandwidth(bw);
	rx_usrp->set_rx_bandwidth(bw);

	tx_usrp->set_tx_freq(freq);
	rx_usrp->set_rx_freq(freq);

	tx_usrp->set_tx_gain(gain);
	rx_usrp->set_rx_gain(gain);
	
	tx_usrp->set_tx_subdev_spec(subdev);
	rx_usrp->set_rx_subdev_spec(subdev);    // tx the same as rx?

	std::cout<<boost::format("Using TX Device: %s") % tx_usrp->get_pp_string()<<endl;
	std::cout<<boost::format("Using RX Device: %s") % rx_usrp->get_pp_string()<<endl;

	// precompute the waveform values
	string wt = "SINE";
	const wave_table_class wave_table(wt,ampl);
	const size_t step = boost::math::iround(wave_freq/tx_usrp->get_tx_rate() * wave_table_len);
	size_t index = 0;


	// channel detect?

	// set transmit streamer
	uhd::stream_args_t stream_args(cpu,wire);
        int num_channels = stream_args.channels.size();
	uhd::tx_streamer::sptr tx0 = tx_usrp->get_tx_stream(stream_args);
   	int spb = tx0->get_max_num_samps()*10;       // why *10?
	//int spb = 15000;
	rbuff.assign(spb,0);
	buff.assign(spb,0);
	cout<<boost::format("spb = %f\n")%spb;
 	
	//vector<complex<float> > buff(spb);
	//vector<complex<float>* > buffs(channel_nums.size(),&buff.front());

	if(channel_nums.size() > 1)
	{
	tx_usrp->set_time_source(pps);

	tx_usrp->set_time_unknown_pps(uhd::time_spec_t(0.0));
	//sleep(1);
	}
	else tx_usrp->set_time_now(0.0);

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

    // set up metadata
	uhd::tx_metadata_t md;
	md.start_of_burst = true;				// set start of burst to true for 1st packet in the chain
	md.end_of_burst = false;                // set end of burst to true for the last packet in the chain
	md.has_time_spec = true;                // set true to send at the time specified by time spec; set false to send immediately
	md.time_spec = uhd::time_spec_t(0.1);   // what's the meaning -- sending time   "tx_usrp->get_time_now() "
	
    // use thread_group to send and receive data at the same time
	string txfile = "tx_out";
    boost::thread_group transmit_thread;
    transmit_thread.create_thread(boost::bind(&transmitter, wave_table, tx0, md,txfile, step, index, num_channels));		// bind: return tranmitter(), no argument needed


// receive begin
#define recv_to_file_args(format) \
	(rx_usrp, format, wire, file, spb, total_num_samps, total_time)
	//recv to file
	if (type == "double") recv_to_file<std::complex<double> >recv_to_file_args("fc64");
	else if (type == "float") recv_to_file<std::complex<float> >recv_to_file_args("fc32");
	else if (type == "short") recv_to_file<std::complex<short> >recv_to_file_args("sc16");
	else throw std::runtime_error("Unknown type " + type);

	// use the global variables to do the cancellation
	int preamble_length = rate/wave_freq;			// for sine wave: a period
	int estimator_length = 20;
	int pilot_length = 400;
	int signal_length = 1e4;
	VectorXf preamble(preamble_length);
	for(int n = 0; n < preamble_length; n ++)
		preamble(n) = wave_table(index += step).real();		  // preamble sequence
	
	
	// seems a little difficult to convert from vector<complex> to VectorXf
	// using Map(rbuff.data(), rbuff.size()) is wrong if rbuff is complex

	VectorXf x(signal_length), y(signal_length);
	for(int i = 0; i < signal_length; i ++)
	{
		x(i) = buff[i].real();
		y(i) = rbuff[i].real();

	}

	//finished
	stop_signal_called = true;
	VectorXf y_clean = dg_sic(x, y, preamble, estimator_length, preamble_length, pilot_length, signal_length, rate);
	std::cout << std::endl << "Done!" << std::endl << std::endl;
	return EXIT_SUCCESS;

}
