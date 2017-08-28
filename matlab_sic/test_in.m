close all
clear all
path = {'/home/flexicon/Downloads/GNURadio/uhd/host/examples/qpsk_transceiver1/build/';
    '/home/flexicon/Downloads/GNURadio/uhd/host/examples/transceiver/build/';
    'D:/ѧϰ/Flexicon/04code/matlab_sic/qpsk_data/'
    };

txname = {'tx_out','tx_out','tx_file'};
rxname = {'qpsk_2M','sine_wave_2M','rx_file'};
id = 3;

fid = fopen(strcat(path{id},txname{id}));              % read tx
a1 = fread(fid,[2,inf],'float');
fclose(fid);

fid = fopen(strcat(path{id},rxname{id}));              % read rx
a2 = fread(fid,[2,inf],'float');
fclose(fid); 

t = 1:length(a1);
figure
plot(t,a1(1,:),':r');
hold on;
t = 1:length(a2);
plot(t,a2(1,:),'b');