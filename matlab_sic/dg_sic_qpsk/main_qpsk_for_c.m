% main_qpsk reading TX and RX from file of qpsk_transceiver1.cpp

close all
clear
path = '';
name = {'tx_qpsk','rx_qpsk'};

%% read TX, RX from file
a = cell(1,2);
for i = 1:2
    fid = fopen(strcat(path,name{i}),'rb');
    a{i} = fread(fid,[2 inf],'float');
    fclose(fid);
end
   
%  1, it is real number and do the cancellation only for real
% tx = a{1}(1,:);     
% rx = a{2}(1,:);

% 2, use the complex number
tx = a{1}(1,:) + 1i*a{1}(2,:);     
rx = a{2}(1,:) + 1i*a{2}(2,:);

figure
plot(1:length(tx),real(tx),':r');
hold on
plot(1:length(rx),real(rx),'b');

%% detect the correct buffer manually
buff_length = 2e4;        % total length of the window, may not all the signal
st1 = 510992;             % st1, st2 only for default file, change with signals
st2 = 517442;
tx_buff = tx(st1:st1 + buff_length - 1);
rx_buff = rx(st2:st2 + buff_length - 1);

%% do the cancellation
rate = 2e6;
sps = 8;
N_T = sps;
signal_length = 1e4;
preamble_length = 128;
estimator_length = 120;
pilot_length = 640;
preamble = tx_buff(1:preamble_length);
pilot = tx_buff(preamble_length + 1:preamble_length + pilot_length);
data_length = signal_length - preamble_length - pilot_length;
start = preamble_length + 1;        % defn: start of pilot, begin estimation

y_clean = dg_sic_qpsk(tx_buff,rx_buff,rate,N_T,preamble,pilot, data_length, estimator_length,start);
