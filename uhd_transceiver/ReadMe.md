### Introduction
These are the C++ codes for UHD and are tested on USRP B205mini. There are two ways to classify them: based on their uses or signals they deal with. 

### Classification by use
#### 1. Transceiver
Transceiver can generate different waveforms, send and receive them with USRP and usually write into file. The first transceiver, transceiver0.cpp, is for sine wave and it comes from a uhd example, txrx_loop_back_to_file.cpp. Generally, transceiver is used to observe the original signal without processing. It can also provide proper TX and RX signals for Matlab code to do the cancellation, like what qpsk\_transceiver1.cpp does.

#### 2. Tranceler
Tranceler or transceiver\_canceler can not only transceive signal but also do the cancellation. There are trancelers for single sine wave, multi-tone sine wave and QPSK signal. And they can be divided into offline canceler and real time canceler (not true "real time" but can update cancellation result continuously). By far, canceler for single sine wave, multi-tone sine wave have been realized basically and tested with USRP. Real time canceler for QPSK signal is realized but doesn't manage to cancel SI due to bad synchronization.

Besides, linear cancellation is tested reliable while all cancelers are updated to non linear cancellation since the latter can be linear by setting dim = 1 defaultly. 

Nonlinear cancellation prove reliable as long as generated sine wave has a precise frequency, which can be only realized based on VectorXf defined signal intead of wavetable. If not, the shift of sample will increase the rank of the matrix A in digital cancellation algorithm, which will influence channel estimation. By now, I have realized VectorXf based funcion of singal generation in tranceler_mt_rt\_final.  

#### 3. Function and header file
digital\_SIC.cpp, rcos\_filter.cpp contain functions. They are now rewritten into every source file.

wavetable.hpp, fft.hpp are header files. The former is to generate different signals and the latter is to realize fft and sic\_db().

sic\_db.cpp is used to calculate the cancellation result in dB through a fft block from github. Thanks for wareya's work. See more at [Public-domain single-header FFT library (power-of-2 size case only)][1]. Now it is in fft.hpp.

nonlinear.cpp is used to generate matrix A with nonlinear components for nonlinear cancellation. It also contains a back up of x2A() for linear cancellation only before.

### Classification by signal
#### 1. single sine wave
transceiver\_0.cpp, transceiver\_canceler.cpp and transceiver_canceler_rt\_sync.cpp are for single sine wave. 

* transceiver\_0.cpp is the first UHD transceiver from UHD example. It sets the USRP, generates sine wave from waveform.hpp and uses buffer to send and receive signal. Receiver, recv_to_file(), is in the main thread while transmitter() is in another thread to make sure they work simultaneously. 

* transceiver\_canceler.cpp has an offline cancellation part in another thread, which does about 50 dB cancellation and doesn't interrupt transmitter or receiver. It writes TX, RX and clean signal into tx\_file, rx\_file and y_clean_file so we can visualize them in Matlab. And it also writes buffers to file in transmitter and receiver function.

* transceiver_canceler_rt\_sync.cpp has a real time cancellation part in another thread. Compared to tranceler above, it do the cancellation without any interruption because it doesn't write any file in transmitter or receiver. It also has a signal synchronization part (though may not work now). 

Besides, digital\_SIC.cpp is the first version of cancellation part in any canceler. It is out of date now, however, because no other code calls it.

#### 2. multi-tone sine wave
transceiver_canceler_multi\_tone.cpp, tranceler_multi_tone\_rt.cpp tranceler_mt_rt\_final.cpp are for multi-tone sine wave.

* transceiver_canceler_multi\_tone.cpp is the first version for multi-tone sine wave, which changes TX from above. It does offline cancellation in another thread for 4 tones now, 100 kHz to 400 kHz with a sampling rate of 2 MHz. It is tested. 

* tranceler_multi_tone\_rt.cpp adds real time and signal synchronization parts and the cancellation hasn't been tested yet. 

* tranceler_mt_rt\_final.cpp, the latest version of cancellation example, has several changes compared to code above: 
  1.  based on VectorXcf instead of vector; 
  2.  move x2A(), greater1() and peaks() to header file nonlinear\_peak.hpp; 
  3.  modify transmitter(), receiver() and digital\_canceler() to make them suitable for multi-tone sine wave;
  4.  generate multi-tone sine wave precisely in function mt_sine_generation().
  
  About running the example, please see the section *How to run the example* below. After make the project, there are two ways to test it. The first is directly typing command
  ```
  ./tranceler_mt_rt_final
  ```
  , and then it will set four-tone sine wave with frequency 100 kHz, 200 kHz, 300 kHz, 400 kHz.
  
  The second is to set the parameters manually. Two formats for setting are provided as following. 
  
  1. Input wave num, wave freq 1 and wave spacing to generate a equidistant multi-tone sine wave, e.g,
  ```
  ./tranceler_mt_rt_final --wave-num 3 --wave-freq-1 200e3 -- wave-space 200e3
  ```
  generates a two-tone sine wave with frequency 200e3, 400e3, 600e3.
  
  2. Input at most 8 tone frequency straightforward.
  ```
  ./tranceler_mt_rt_final --freq-1 200e3 --freq-2 400e3 --freq-3 600e3
  ```
  generates the same multi-tone sine wave as 1).
  
  The program displays information about canceler setting (sampling rate, start sample, length of estimator, pilot, preamble, signal and TX/RX buffer) as well as continuous update on cancellation result. Every certain number of canceled buffers, it shows strength of the specific tone of received signal, clean signal and their difference (cancellation result) in dB. It also calculates average cancellation for the whole signal.


#### 3. QPSK signal
qpsk_tranceler.cpp, qpsk_tranceler_rt.cpp, qpsk_transceiver\_1.cpp and rcos\_filter.cpp are for QPSK signal.

* qpsk\_tranceler.cpp fails to do the cancellation for QPSK signal because it has no synchronization now. 

* qpsk_tranceler_rt.cpp integrates real time part and calculation of cancellation (sic\_db()) but by now it cannot do the right synchronization. 

* qpsk_transceiver_1.cpp can generate and transceive QPSK signal so that we can use them in Matlab and do the cancellation. It calls wave\_generation() and rcos\_filter() (inside qpsk_transceiver_1.cpp itself but from rcos\_filter.cpp) to generate continuous QPSK waveform from complex bipolar code like 1 + 1i, -1 - 1i. It works fine with Matlab cancellation and can get about 30 dB cancellation.

* rcos\_filter.cpp contains primary function about waveform generation of band-limited signal, such as rcos_filter() and wave_generation(). The former generates a raised cosine filter and the latter uses it to generate actual TX waveform like QPSK TX signal. main() generates QPSK TX signal so it can be easily integrated into transceiver and canceler code.


### How to run the example
These examples are all based on UHD TX and RX examples and it is simple to run if you have experience with UHD. And it is highly recommended to take a look at the UHD tutorials. Take linux for example. First, check dependents of the code, Eigen files, boost and uhd library. Then, put these file into the same folder of UHD examples. When you want to run an example, create a new folder, put source, header file and CMakeLists.txt into the folder and then build the project.

There's an folder for QPSK tranceler. It contains fft.hpp, wavetable.hpp and can be build, make and run by these commands:

```
mkdir build
cd build
cmake ..
make
./qpsk_tranceler_rt
```

Notice that if your library or Eigen need root to access, you need to use "sudo make" instead of "make" to make the program.

[1]:https://github.com/wareya/fft
