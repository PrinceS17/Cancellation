### Introduction
This is folder of Matlab code of digital SI cancellation for the full-duplex system. I use Matlab to not only simulate SIC to verify the digital cancellation algorithm but also visualize and test the results from USRP. Basically, they can be classified into 3 parts: SIC for sine wave, SIC for QPSK signal and test code for TX, RX from USRP. 

For the detail of digital SIC algorithm, please refer to the digital cancellation part of [*Full Duplex Radio*][1]\[1\].

### SIC for sine wave
This part contains code for the simulation and test for single sine wave and is located in *dg_sic_sine* folder. Main source files are dgSIC_sine_0.m, dgSIC_sine_rx\_1.m and dg\_sic.m, mat\_generation.m are functions that most related to them.

* digSIC_sine_0.m is for the simulation for single sine wave. It generates single sine wave as TX and sine wave with AWGN and ISI as RX. By running it, we can test the algorithm and get about 50 dB linear cancellation.

* dig_sine_rx\_1.m is for the test of SIC algorithm in Matlab. It obtains TX and RX signals from UHD transceiver for single sine wave and then calls dg\_sic() and mat\_generation() to do the cancellation. It can both do the linear and the nonlinear cancellation by changing the parameter *dim*. Since it is for offline cancellation, the program only works for TX/RX signal with length of *signal_length*. A signal consists of preamble(not necessary for sine wave), pilot and data and I use an estimator which scans from i-k to i+k-1. Their length are defined as *estimator_length, pilot_length*.   

### SIC for QPSK signal


### Test code for TX, RX from USRP
Source files that have prefix "test" are test code for signals from USRP like TX, RX and cleaned signal. Generally, they read from a path where UHD program writes data into file and then visualize the cancellation result or TX/RX signals. 

* test\_in.m is only for visualization of initial TX and RX signal. Since it is harder to find the problem after the program processes the signal, I use this code to observe the original signal. However, if you want to visualize a signal in your own path, you have to modify the path. What I usually do is to add an option to the path, txname, rxname arrays. 

  Notice that it reads data in 2 rows for most UHD programs write complex signal into file before. And it only plots the real part of the signal by using 1st row of a1 and a2, which can be easily changed.

* test_trans_canc.m is for visualization of pilot, data and spectrum of TX, RX and cleaned signal (x, y, y\_clean). Notice that it reads data in 1 row since UHD cancelers write only real data into file. Like above, you need to modify the build\_path and filename to add new data file to visualize.

* test\_trans.m is early version of test_trans_canc.m, which can only plot data and spectrum of sine wave.

* test_trans_canc\_rt.m is for continuous visualization of spectrum of TX, RX and cleaned signal and I often use it to see the result of real time cancellation. Settings for build\_path and filename are the same as above.

### Other functions
* figure\_sic.m is used to plot all the TX and RX results: channel h, preamble, pilot, data and spectrum. It is often called by simulation code where preamble and h can be easily obtained.

### Reference
[1] Bharadia, D., McMilin, E., & Katti, S. (2013). Full duplex radios. ACM SIGCOMM Computer Communication Review, 43(4), 375-386.

[1]:https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwjI4K7U1JHWAhVH9IMKHR2xD0sQFghuMAA&url=https%3A%2F%2Fwww.stanford.edu%2F~skatti%2Fpubs%2Fsigcomm13-fullduplex.pdf&usg=AFQjCNGZDqwpXhxTrJmdkXovcJt1N28TkQ

