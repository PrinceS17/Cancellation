% received signal from sine wave, do the cancellation
clear
close all

%% parameter definition
rate = 2e6;
freq = 100e3;
signal_length = 1e4;
estimator_length = 42;
estimator_start = 30;
pilot_length = 999;
data_length = 200;
t = (1:signal_length)/rate;
dim = 1;

%% obtain the transmitted signal also from file
% x = 1*sin(2*pi*freq*t);
% 
%  fid = fopen('tx_out');
% a = fread(fid, [2, inf], 'float');
% x = a(1,1:signal_length) + 1i*a(2,1:signal_length);
% x = real(x);
% fclose(fid);
% 
% %% obtain the received SINE wave from data file
% fid = fopen('sine_wave_M');
% a = fread(fid, [2, inf], 'float');
% st = 2e6;
% y = a(1,st:st + signal_length - 1) + 1i*a(2,st:st + signal_length - 1);
% y = real(y);     % only use real part now
% fclose(fid);

load xyz.mat
x = pilot_I';
y = v_pilot_avg_d;

%% digital cancellation
y_clean = dg_sic(x,y,rate,freq,data_length,estimator_length,estimator_start,pilot_length,signal_length,dim);

%% plot: time & frequency domain
fc = 400e3;
bw = 10;
rg = 5e3;
sic_db([y(1:length(y_clean));y_clean], rate, fc, bw, rg),
% plot_sic(x,y,y_clean,rate);

