% received signal from sine wave, do the cancellation

%% parameter definition
rate = 2e6;
freq = 100e3;
signal_length = 1e4;
estimator_length = 80;
estimator_start = 50;
pilot_length = 1000;
t = (1:signal_length)/rate;

%% obtain the transmitted signal also from file
 x = 1*sin(2*pi*freq*t);

%  fid = fopen('tx_out');
% a = fread(fid,[2 inf], 'float');
% x = a(1,1:signal_length) + 1i*a(2,1:signal_length);
% x = real(x);
% fclose(fid);

%% obtain the received SINE wave from data file
fid = fopen('sine_wave_2M');
a = fread(fid, [2, inf], 'float');
st = 1;
y = a(1,st:signal_length) + 1i*a(2,st:signal_length);
y = real(y);     % only use real part now
fclose(fid);

%% digital cancellation
y_clean = dg_sic(x,y,rate,freq,estimator_length,estimator_start,pilot_length,signal_length);

%% plot: time & frequency domain
% plot_sic(x,y,y_clean,rate);


