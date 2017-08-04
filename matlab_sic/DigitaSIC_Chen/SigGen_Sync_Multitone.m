clc;
clear;
close all;

pilot = 1/10; % fraction of the signal that is the pilot

load('LabViewData_Sync_Multitone.mat') % load the raw LabView data
% T_I_Rx = time values for received in phase waveform
% Sig_I_Rx = received in phase waveform
% T_Q_Rx = time values for received quadrature waveform
% Sig_Q_Rx = received quadrature waveform
% T_I_Tx = time values for transmitted in phase waveform
% Sig_I_Tx = transmitted in phase waveform
% T_Q_Tx = time values for transmitted quadrature waveform
% Sig_Q_Tx = transmitted quadrature waveform


pilot_I = Sig_I_Tx(1:length(Sig_I_Tx)*pilot);
pilot_Q = Sig_Q_Tx(1:length(Sig_Q_Tx)*pilot);
sig_I = Sig_I_Tx(length(Sig_I_Tx)*pilot+1:length(Sig_I_Tx));
sig_Q = Sig_Q_Tx(length(Sig_Q_Tx)*pilot+1:length(Sig_Q_Tx));
N_tot_pilot = pilot*length(Sig_I_Tx);
N_tot = length(Sig_I_Tx);
N_tot_signal = (1-pilot)*length(Sig_I_Tx);

% Take the FFT of the pilot and the received signal, multiply, then take
% the ifft in order to get the pilot convolved with the received signal.
% fft_pilot = fft(pilot_I,10000);
% fft_rx = fft(Sig_I_Rx,10000);
% ifftconv = ifft(fft_pilot.*fft_rx,10000);
% plot(ifftconv);
% figure;
% [maxval, index] = max(ifftconv);
% index = index + 10;

sig_I_conv = conv(pilot_I, Sig_I_Rx)
[maxval, index] = max(sig_I_conv)
figure(10); hold on;

% Use the index to shift the received signal the right amount (the received
% signal is renamed vout)
vout = [transpose(Sig_I_Rx(index-length(pilot_I)+1:length(Sig_I_Rx))) transpose(Sig_I_Rx(1:index-length(pilot_I)))];
plot(vout);
rx_pilot_I = vout(1:N_tot_pilot)

