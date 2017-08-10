close all
clear all

%% definition
samp_rate = 3e6;
estimator_length = 60;
pilot_length = 400;
preamble_length = 64;
signal_length = 5e3;
samp_per_sym = 8;
ampl = 0.3;
gain = 20;
t = (1:signal_length)/samp_rate;
start = preamble_length + estimator_length/2 + 1;

%% wave generation
symbol = [1 + 1i, 1 - 1i, -1 + 1i, -1 - 1i];
preamble_sym = (1 + 1i)*ones(1,preamble_length/samp_per_sym);          % generate symbol sequence
pilot_sym = symbol(randi(4,1,pilot_length/samp_per_sym));

tx_data = symbol(randi(4,1,(signal_length - preamble_length - pilot_length)/samp_per_sym));
x = [preamble_sym, pilot_sym, tx_data];
x = qpsk_generation(x,0.5,4,samp_per_sym);

%% write tx data & read received data from USRP
% write tx data into file
% addpath(genpath('/home/flexicon/Downloads/GNURadio/uhd/host/examples/qpsk_transceiver/build'));
% fid = fopen('qpsk_tx','wb');
% xw = reshape([real(x);imag(x)],1,length(x)*2);
% fwrite(fid, xw,'float');
% fclose(fid);

% % 1, generate received signal for simulation
% phi = pi/7;
% delay = 100;
% SNR = 50;
% y = [zeros(1,delay), x*exp(1i*phi)];             % y is longer than x
% y = awgn(y,SNR,'measured');                      % add noise

% 2, use received data from USRP
fid = fopen('qpsk_rx_file','rb');
a = fread(fid, [2, inf], 'float');               
% tx_beg = 2.407e5;
tx_beg = 1;
tx_end = tx_beg + signal_length - 1;
y = a(1, tx_beg:tx_end) + 1i*a(2,tx_beg:tx_end);
fclose(fid);

%% digital cancellation and figures
preamble = x(1:preamble_length);
pilot = x(preamble_length + 1:preamble_length + pilot_length);
y_clean = dg_sic_qpsk(x,y,samp_rate,samp_per_sym,preamble,pilot,estimator_length,start);
