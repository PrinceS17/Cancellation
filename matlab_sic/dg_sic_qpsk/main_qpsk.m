% digital cancellation with QPSK modulation
% 1, how to design preamble and pilot
% 2, how to set the critieria for cancellation

%% parameter definition & wave generation
rate = 1e6;                                     % symbol rate
estimator_length = 20;                          % these number are all of symbol 
pilot_length = 400;
preamble_length = 32;                           % a preamble of 4 byte
signal_length = 2e4;
samp_per_sym = 8;
ampl = 0.3;
gain = 20;
t = (1:signal_length*samp_per_sym)/rate;
start = preamble_length + estimator_length/2 + 1;

symbol = [1 + 1i, 1- 1i, -1 + 1i, -1 - 1i];
preamble = (1 + 1i)*ones(1,preamble_length);                                    % for synchronization
% pilot = reshape(symbol'*ones(1,round(pilot_length/4)),1,pilot_length);          % high frequency: not work
% pilot = reshape(ones(round(pilot_length/4),1)*symbol,1,pilot_length);           % low frequency: worked, MSE = 2e-28, bias 1e-14
pilot = symbol(randi(4,1,pilot_length));                                      % all frequency: random, work, MSE = 2e-30, bias 1e-15

tx_data = symbol(randi(4,1,signal_length - preamble_length - pilot_length));    % TX data we know
x = [preamble, pilot, tx_data];                                                 % transmit signal
x = qpsk_generation(x,0.5,4,samp_per_sym);


%% write tx data & read received data from USRP
% fid = fopen('qpsk_output','wb');
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
fid = fopen('qpsk_tx','rb');
a = fread(fid, [2, inf], 'float');               
tx_beg = 2.407e5;
tx_end = tx_beg + signal_length*samp_per_sym - 1;
y = a(1, tx_beg:tx_end) + 1i*a(2,tx_beg:tx_end);
fclose(fid);
% figure(5);plot(real(y));


%% digital cancellation
bound1 = preamble_length*samp_per_sym;
bound2 = bound1 + pilot_length*samp_per_sym;
preamble_seq = x(1:bound1);
pilot = x(bound1 + 1:bound2);
estimator_length_seq = estimator_length*samp_per_sym;
start = start*8 - 7;
[y_clean, MSE] = dg_sic_qpsk(x, y, rate*samp_per_sym, preamble_seq, pilot, estimator_length_seq, start);
MSE,

%% plot: time, frequency domain and constellation
figure(1);
plot_sic(real(x), real(y), real(y_clean), t, start, signal_length, 6*rate);

figure(2);
plot_constellation(x, y, y_clean);

