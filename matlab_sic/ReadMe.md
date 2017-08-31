### Introduction
This is folder of Matlab code of digital SI cancellation for the full-duplex system. I use Matlab to not only simulate SIC to verify the digital cancellation algorithm but also visualize and test the results from USRP. Basically, they can be classified into 3 parts: SIC for sine wave, SIC for QPSK signal and test code for TX, RX from USRP.

#### SIC for sine wave


#### SIC for QPSK signal


#### Test code for TX, RX from USRP
Source files that have prefix "test" are test code for signals from USRP like TX, RX and cleaned signal. Generally, they read from a path where UHD program writes data into file and then visualize the cancellation result or TX/RX signals. 
