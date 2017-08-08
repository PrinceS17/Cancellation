% test tranceiver.cpp
close all
clear all

addpath(genpath('/home/flexicon/Downloads/GNURadio/uhd/host/examples/transceiver_canceler/build'));
fid = fopen('tx_file');               % read x
x = fread(fid,[1,inf],'float');
fclose(fid);

fid = fopen('rx_file');            % read y
y = fread(fid,[1,inf],'float');
fclose(fid);

fid = fopen('y_clean_file');              % read y_clean
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
ylim([-0.25 0.25]); xlim([0 1e-4]);

plot_sic(x,y,y_clean,rate);
