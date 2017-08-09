% test tranceiver_canceler
close all
clear all
build_path = {'/home/flexicon/Downloads/GNURadio/uhd/host/examples/transceiver_canceler/build';
    '/home/flexicon/Downloads/GNURadio/uhd/host/examples/transceiver_canceler_multi_tone/build'};
filename = {'tx_file','rx_file','y_clean_file';
    'tx_file_mt','rx_file_mt','y_clean_file_mt'};

id = 1;           % switch: 1 for transceiver_canceler; 2 for multi_tone

addpath(genpath(build_path{id}));
fid = fopen(filename{id,1});              % read x
x = fread(fid,[1,inf],'float');
fclose(fid);

fid = fopen(filename{id,2});              % read y
y = fread(fid,[1,inf],'float');
fclose(fid);

fid = fopen(filename{id,3});              % read y_clean
y_clean = fread(fid,[1,inf],'float');
fclose(fid);
% subplot(2,1,1);
% plot(a1(1,:),':r'); hold on
% plot(a2(1,:));
% legend('y','y_{clean}');

% plot
rate = 2e6;   % should consistent with the C++ code!!!
f = 100e3;
L = length(y);
t = (1:L)/rate;
% x = 0.5*sin(2*pi*f*t);
figure
plot(t,x,':r',t,y,'b'); hold on;
plot(t(1:length(y_clean)),y_clean,'g'); % xlim([0 8e-5]); ylim([-0.1 0.1]);
% aksdjflkajf;jsaylim([-0.25 0.25]); xlim([0 1e-4]);

plot_sic(x,y,y_clean,rate);
