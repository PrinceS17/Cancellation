transceiver_0.cpp is now mainly for SINE wave and it uses wavetable.hpp to generate all kinds of waveforms. 

digital_SIC.cpp is about the cancellation part which comes from the following part.

transceiver_canceler.cpp is the combination of transceiver and canceler and is tested with USRP. It can do the cancellation for sine wave by reading a buff after transmitting and receiving finish and write the x, y, y_clean to tx_file, rx_file and y_clean_file. So now it's the offline cancellation. Then we can use Matlab to read the signal and visualize the cancellation result.
