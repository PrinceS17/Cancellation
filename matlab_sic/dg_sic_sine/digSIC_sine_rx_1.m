% received signal from sine wave, do the cancellation

%% parameter definition
rate = 2e6;
freq = 100e3;
signal_length = 1e4;
estimator_length = 80;
estimator_start = 50;
pilot_length = 1000;
t = (1:signal_length)/rate;
x = 1*sin(2*pi*freq*t);

%% obtain the received SINE wave from data file
fid = fopen('sine_wave_2M');
a = fread(fid, [2, inf], 'float');
y = a(1,:) + 1j*a(2,:);
y = y(1:signal_length);
y = real(y);                    % only use real part now

%% digital cancellation
y_clean = dg_sic(x,y,estimator_length,estimator_start,pilot_length,signal_length);

%% plot: time & frequency domain
plot_sic(x,y,y_clean,t,estimator_start,signal_length,rate);


