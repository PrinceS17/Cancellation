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
#include <time.h>
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

// real time cancellation for single sine wave


static bool stop_signal_called = false;
int global_num[2] = {0,0};			 // num of buffer in global_rx, tx
vector<complex<float> > global_rx;		 // global_rx, tx has many buffers, like 1000
vector<complex<float> > global_tx;	
int gsize;
// int tx_pos = 0;					 // position: index of latest new entry of global_tx
// int rx_pos = 0;					 // like above
boost::mutex mtx;					 // mtx for reading global variables

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

int dg_sync(VectorXf& preamble, VectorXf& rbuff)   // return delay in rbuff
{
	int cor_length = rbuff.size() - preamble.size() + 1;
	if(cor_length < 0)
	{
		cout<<"-- error: cor length is negative!"<<endl;
		exit(0);
	}
	VectorXf Cor(cor_length);					   // make the preamble
	
	
	for(int i = 0; i < cor_length; i++) 
		Cor(i) = preamble.transpose()*rbuff.segment(i,preamble.size());
	/*Index idx;					        // simply find the max number
	Cor.maxCoeff(&idx);*/
	int peak_num = 7;
	VectorXf idxs = peaks(Cor, peak_num);				// find the biggest 10 peaks
	sort(idxs.data(),idxs.data() + peak_num,less<int>());		// find the earliest idx among the 10
	return idxs(0);
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

	// write estimated pilot into file	
	ofstream out;
	out.open("estimated_pilot",ios::out | ios::binary);
	VectorXf yp = A*h;
	out.write((const char*)yp.data(),yp.size()*sizeof(float));
	out.close();

	out.open("rx_pilot",ios::out | ios::binary);
	out.write((const char*)rbuff.data(),rbuff.size()*sizeof(float));
	out.close();
	
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
	//VectorXf &x, 			// no need when using global variables
	//VectorXf &y,			              // initial signal got from UHD: here haven't defined complex number
	int spb,	
	VectorXf &preamble_vf,
	VectorXf &preamble,		              // should have complete definition later
	int estimator_length,
	int preamble_length,
	int pilot_length,
	int signal_length,			      // the length of TX which made signal_length + delay <= RX's length, may be redefine
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
	
	ofstream outfile[3];
	string name[3] = {"tx_file","rx_file","y_clean_file"};
	int num[3] = {0,0,0};
	int can_num = 0;
	int delay_const = -1;

	for(int i = 0; i < 3; i++)
		outfile[i].open(name[i].c_str(), ios::out | ios::binary);
	
	VectorXf x(signal_length), y(signal_length);
	VectorXf rx1;
	bool if_synced = 0;	

	while(not stop_signal_called)
	{
	

		mtx.lock();	// begin reading global buffer
		int tx_num = global_num[0];
		int rx_num = global_num[1];
		if(!if_synced && rx_num > 3)		// for sync: suppose 5 buffer is the biggest delay
		{
			cout<<"-- synced! "<<endl;
			rx1 = VectorXf::Zero(global_num[1]*spb);
			for(int i = 0; i < rx1.size(); i ++)
				rx1[i] = global_rx[i % gsize].real();
			if_synced = 1;
		}
		else if(!if_synced) { mtx.unlock(); continue;}			// before synchronization
		else ;

		if(delay_const > 0) 						// begin cancellation after 1st sync
		{
			int can_pos = can_num*signal_length;
				
			for(int i = 0; i < signal_length; i ++)
			{
				x(i) = global_tx[(can_pos + i) % gsize].real();
				y(i) = global_rx[(can_pos + i + delay_const) % gsize].real();	// detail of y_clean's length should be notice
			}
			can_num ++;							// record the number of cancellation
		}
		mtx.unlock();	// end reading global buffer

				
		if(delay_const < 0)						// if rx1 read but not synchronized
		{
			delay_const = dg_sync(preamble_vf, rx1);		// synchronize for TX and RX			
			cout<<"-- TX, RX delay = "<<delay_const<<endl;	
			continue;
 		}
		if(can_num%100 == 0)
			cout<<"-- TX No. "<<tx_num<<" , RX No. "<<rx_num<<" , Cancel No. "<<can_num*signal_length/spb<<endl;
	
		int k = estimator_length/2;
		int delay = dg_sync(preamble, y);
		//cout<<"-- delay = "<<delay<<endl;
	
		// define tx&rx_pilot and estimate h
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
		VectorXf y_clean = dg_cancel(tx_data, rx_data, h, estimator_length);

		float * ptr[3] = {x.data(), y.data(), y_clean.data()};
		int length[3] = {x.size(), y.size(), y_clean.size()};

		// write to file and some other work
		for (int i = 0; i < 3; i++)
		{
			outfile[i].write((const char*)ptr[i], length[i]*sizeof(float));     // try I/Q channel by one first
			num[i] += length[i];

			if(num[i]*sizeof(float) > 50e6*sizeof(char)) 
			{					
				outfile[i].close();
				outfile[i].open(name[i].c_str(),ios::binary | ios::out);
				num[i] = 0;
				cout<<"-- New "<<name[i]<<endl;
			}
		}
	}
	for(int i = 0; i < 3; i++)
		outfile[i].close();
}



// tranceiver begin

void transmitter(
    //std::vector <complex<float> > &buff,	// no need, transceiver can also generate buff inside the function
    int spb,
    wave_table_class wave_table,
    vector<complex<float> > preamble_1,
    uhd::tx_streamer::sptr tx0,
    uhd::tx_metadata_t md,
const string file,
    size_t step[4],
    size_t index[4],
    int num_channels
)
{

	int num = 0;
	vector <complex<float> > buff(spb), pre_buff(spb);	// buffer: sending to USRP, pre: calculate before	
	
	// calculate pre-buff outside the loop to save resources
	
	for(size_t n = 0; n < buff.size(); n ++)
	{	
		pre_buff[n] = (float)0.25*wave_table(index[0]) + (float)0.25*wave_table(index[1])
			+ (float)0.25*wave_table(index[2]) + (float)0.25*wave_table(index[3]);
		for(int i = 0; i < 4; i ++)	index[i] += step[i];
	}

	int if_preamble = 1;
	while(not stop_signal_called)	
	{
		if(if_preamble)
		{
			buff = preamble_1;
			if_preamble --;
			cout<<"-- if preamble? "<<if_preamble<<endl;
		}
		else buff = pre_buff;

		mtx.lock();
		int tx_pos = spb*global_num[0];
		for(int i = 0; i < spb; i ++)
			global_tx[(tx_pos + i) % gsize] = buff[i];
		
		global_num[0] ++;
		//cout<<"-- tx num: "<<global_num[0]<<endl;	
		mtx.unlock();

		num += tx0->send(&buff.front(),buff.size(),md);   // send tx data from buff
	
		md.start_of_burst = false;
		md.has_time_spec = false;

	}
	
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
	std::vector<samp_type> rbuff(samps_per_buff);    
	
	//string output;
	unsigned long long num_total_samps = 0;
	int spb = samps_per_buff;

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
	
		mtx.lock();
		int rx_pos = spb*global_num[1];
		for(int i = 0; i < spb; i ++)
			global_rx[(rx_pos + i) % gsize] = rbuff[i];
		global_num[1] ++;
		//cout<<"-- rx num: "<<global_num[1]<<endl;
		mtx.unlock();

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

	}

	// end receiveing

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

	double rate,rx_rate,tx_rate,freq,gain,bw,ampl;
	rate = 2e6;
		file = file +"_" + boost::lexical_cast<string>(rate/1e6) + "M";
	double wave_freq[4] = {100e3, 200e3, 300e3, 400e3};
	freq = 915e6;
	gain = 20;				// loop cable: 25; w/o cable: 45
	bw = 1e6;
	rx_rate = rate;
	tx_rate = rate;
	ampl = 0.2;
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
	string pre_wt = "SINE";
	const wave_table_class wave_table(wt,ampl);
	const wave_table_class preamble_table(pre_wt,2*ampl);
	size_t step[4];
	for(int i = 0; i < 4; i ++)	
		 step[i] = boost::math::iround(wave_freq[i]/tx_usrp->get_tx_rate() * wave_table_len);
	size_t index[4] = {0,0,0,0};

	
	// channel detect?

	// set transmit streamer
	uhd::stream_args_t stream_args(cpu,wire);
        int num_channels = stream_args.channels.size();
	uhd::tx_streamer::sptr tx0 = tx_usrp->get_tx_stream(stream_args);
   	int spb = tx0->get_max_num_samps()*10;       // why *10?
	//int spb = 15000;
	
	int num_buff = 1000;				// test which is fine enough
	global_rx.assign(num_buff*spb,0);		// initialize the global variables
	global_tx.assign(num_buff*spb,0);
	gsize = global_tx.size();
	cout<<boost::format("spb = %f\n")%spb;
 	
	// generate 1st preamble for TX, RX synchronization
	vector <complex<float> > preamble_1(spb);	
	int id_temp = 0;
	for(int i = 0; i < spb; i ++)
		preamble_1[i] = preamble_table(id_temp += step[0]);	// use 100kHz as 1st preamble

	
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
    boost::thread_group threads;
	
    threads.create_thread(boost::bind(&transmitter, spb, wave_table, preamble_1, tx0, md,txfile, step, index, num_channels));		// bind: return tranmitter(), no argument needed

// use the global variables to do the cancellation
	int preamble_length = rate/wave_freq[0];			// for sine wave: a period
	int estimator_length = 60;
	int pilot_length = 400;
	int signal_length = spb;				// control the buffer of a cancellation
	
	VectorXf preamble(preamble_length);
	for(int i = 0; i < 4; i++)
		index[i] = 0;
	for(int n = 0; n < preamble_length; n ++)
	for(int i = 0; i < 4; i ++)
		preamble(n) = wave_table(index[i] += step[i]).real();		  // preamble sequence
	
	VectorXf preamble_vf(2*spb);			// 1st preamble for dg_sic()
	for(int i = 0; i < 2*spb; i ++)
		preamble_vf[i] = preamble_1[i%spb].real();
	cout<<"-- after assign..."<<endl;

	
	threads.create_thread(boost::bind(&dg_sic, spb, preamble_vf, preamble, estimator_length, preamble_length, pilot_length, signal_length, rate));


// receive begin
#define recv_to_file_args(format) \
	(rx_usrp, format, wire, file, spb, total_num_samps, total_time)
	//recv to file: the 1st line required global_buff = buff to be complex<double> which is complicated for now

	/*if (type == "double") recv_to_file<std::complex<double> >recv_to_file_args("fc64");
	else if (type == "float") recv_to_file<std::complex<float> >recv_to_file_args("fc32");
	else if (type == "short") recv_to_file<std::complex<short> >recv_to_file_args("sc16");
	else throw std::runtime_error("Unknown type " + type);*/

	// because now recv is not a thread but in the main, so it has to be the end
		
	recv_to_file<std::complex<float> >recv_to_file_args("fc32");	
	
	
	// seems a little difficult to convert from vector<complex> to VectorXf
	// using Map(rbuff.data(), rbuff.size()) is wrong if rbuff is complex

	/*VectorXf x(signal_length), y(signal_length);
	for(int i = 0; i < signal_length; i ++)
	{
		x(i) = buff[i].real();
		y(i) = rbuff[i].real();

	}*/ // only for cancellation offline

	//finished
	stop_signal_called = true;	
	std::cout << std::endl << "Done!" << std::endl << std::endl;
	return EXIT_SUCCESS;

}
