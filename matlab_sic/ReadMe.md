### Introduction
This is folder of Matlab code of digital SI cancellation for the full-duplex system. I use Matlab to not only simulate SIC to verify the digital cancellation algorithm but also visualize and test the results from USRP. Basically, they can be classified into 3 parts: SIC for sine wave, SIC for QPSK signal and test code for TX, RX from USRP.

#### SIC for sine wave


#### SIC for QPSK signal


#### Test code for TX, RX from USRP
Source files that have prefix "test" are test code for signals from USRP like TX, RX and cleaned signal. Generally, they read from a path where UHD program writes data into file and then visualize the cancellation result or TX/RX signals. 

* test\_in.m is only for visualization of initial TX and RX signal. Since it is harder to find the problem after the program processes the signal, I use this code to observe the original signal. However, if you want to visualize a signal in your own path, you have to modify the path. What I usually do is to add an option to the path, txname, rxname arrays. 

Notice that it reads data in 2 rows for most UHD programs write complex signal into file before. And it only plots the real part of the signal by using 1st row of a1 and a2, which can be easily changed.

* test_trans_canc.m is for visualization of pilot, data and spectrum of TX, RX and cleaned signal (x, y, y\_clean). Notice that it reads data in 1 row since UHD cancelers write only real data into file. Like above, you need to modify the build\_path and filename to add new data file to visualize.

* test\_trans.m is early version of test_trans_canc.m, which can only plot data and spectrum of sine wave.

* test_trans_canc\_rt.m is for continuous visualization of spectrum of TX, RX and cleaned signal and I often use it to see the result of real time cancellation. Settings for build\_path and filename are the same as above.

#### Other functions
* figure\_sic.m is used to plot all the TX and RX results: channel h, preamble, pilot, data and spectrum. It is often called by simulation code where preamble and h can be easily obtained.


