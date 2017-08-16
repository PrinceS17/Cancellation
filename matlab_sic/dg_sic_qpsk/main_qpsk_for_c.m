% main_qpsk reading TX and RX from file of qpsk_transceiver1.cpp

close all
clear all
path = '/home/flexicon/Downloads/GNURadio/uhd/host/examples/qpsk_transceiver1/build/';
name = {'tx_out','qpsk_2M'};

%% read TX, RX from file
a = cell(1,2);
for i = 1:2
    fid = fopen(strcat(path,name{i}),'rb');
    a{i} = fread(fid,[2 inf],'float');
    fclose(fid);
end
    
tx = a{1}(1,:);      % it is real number and do the cancellation only for real
rx = a{2}(1,:);

% figure
% plot(1:length(tx),tx,':r');
% hold on
% plot(1:length(rx),rx,'b');

%% detect the correct buffer manually
buff_length = 1.5e4;        % total length of the window, may not all the signal
st1 = 735842;
st2 = 735842;
tx_buff = tx(st1:st1 + buff_length - 1);
rx_buff = rx(st2:st2 + buff_length - 1);

%% do the cancellation
rate = 2e6;
sps = 8;
N_T = sps;
signal_length = 1e4;
preamble_length = 128;
estimator_length = 80;
pilot_length = 640;
preamble = tx_buff(1:preamble_length);
pilot = tx_buff(preamble_length + 1:preamble_length + pilot_length);
data_length = signal_length - preamble_length - pilot_length;
start = preamble_length + 1;        % defn: start of pilot, begin estimation

y_clean = dg_sic_qpsk(tx_buff,rx_buff,rate,N_T,preamble,pilot, data_length, estimator_length,start);