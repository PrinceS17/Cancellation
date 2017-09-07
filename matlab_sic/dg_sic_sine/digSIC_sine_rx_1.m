% received signal from sine wave, do the cancellation
clear
close all

%% parameter definition
rate = 2e6;
freq = 100e3;
signal_length = 1e4;

pilot_length = 400;
estimator_length = 26;
estimator_start = 38;

data_length = 7000;
t = (1:signal_length)/rate;
dim = 1;

%% obtain the transmitted signal also from file

load sine_2M.mat
start = 5e5;
x = tx_sine(1, start + 1:start + signal_length);
y = rx_sine(1, start + 1:start + signal_length);

%% digital cancellation
y_clean = dg_sic(x,y,rate,freq,data_length,estimator_length,estimator_start,pilot_length,signal_length,dim);

%% plot: time & frequency domain
fc = 400e3;
bw = 10;
rg = 5e3;
sic_db([y(1:length(y_clean));y_clean], rate, fc, bw, rg),
% plot_sic(x,y,y_clean,rate);

