### Introduction
These are the C++ codes for UHD and are tested on USRP B205mini. There are two ways to classify them: based on their functions or signals they deal with. 

### Classification by functions
#### 1. Transceiver
Transceiver can generate different waveforms, send and receive them with USRP and usually write into file. The first transceiver, transceiver0.cpp, is for sine wave and it comes from a uhd example, txrx\_loop\_back\_to\_file.cpp. Generally, 


####　transceiver\_0.cpp 
It is mainly for SINE wave and it uses wavetable.hpp to generate all kinds of waveforms. 

#### digital_SIC.cpp 
It is about the cancellation part which comes from the following part.

transceiver_canceler.cpp is the combination of transceiver and canceler and is tested with USRP. It can do the cancellation for sine wave by reading a buffer after transmitting and receiving finish and write the x, y, y_clean to tx_file, rx_file and y_clean_file. Then we can use Matlab to read the signal and visualize the cancellation result. Note that now it's nearly real time cancellation but without precise synchronization.

transceiver_canceler_multi_tone.cpp is the version for multi-tone sine wave, which changes TX from above.

rcos_filter.cpp contains primary function about waveform generation of band-limited signal, such as rcos_filter() and wave_generation(). The former generates a raised cosine filter and the latter uses it to generate actual TX waveform like QPSK TX signal. main() generates QPSK TX signal so it can be easily integrated into transceiver and canceler code.
