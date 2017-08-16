% test tranceiver_canceler
close all
clear all
addpath(genpath('/home/flexicon/Documents/Cancellation/Cancellation'));
build_path = {'/home/flexicon/Downloads/GNURadio/uhd/host/examples/transceiver_canceler/build';
    '/home/flexicon/Downloads/GNURadio/uhd/host/examples/transceiver_canceler_multi_tone/build';
    '/home/flexicon/Downloads/GNURadio/uhd/host/examples/qpsk_tranceler/build'
    };
filename = {'tx_file','rx_file','y_clean_file','estimated_pilot','rx_pilot';
    'tx_file_mt','rx_file_mt','y_clean_file_mt','estimated_pilot_mt','rx_pilot_mt';
    'tx_file','rx_file','y_clean_file','estimated_pilot','rx_pilot'
    };

id = 3;           % switch: 1 for transceiver_canceler; 2 for multi_tone; 3 for qpsk tranceler

[~,w] = size(filename);
a = cell(1,w);
addpath(genpath(build_path{id}));

for k = 1:w
fid = fopen(filename{id,k});              % read x, y, y_clean and pilot
a{k} = fread(fid,[1,inf],'float');
fclose(fid);
end

x = a{1};
y = a{2};
y_clean = a{3};
est_pilot = a{4};
rx_pilot = a{5};

% plot pilot
rate = 2e6;   % should consistent with the C++ code!!!
f = 100e3;
L = length(y);
t = (1:L)/rate;

if id == 1
    preamble_length = 20;
    pilot_length = 400;
elseif id == 2
    preamble_length = 20;
    pilot_length = 400;
elseif id == 3
    preamble_length = 128;
    pilot_length = 640;
end

figure
id_pilot = preamble_length + 1:preamble_length + pilot_length; 
plot(t(id_pilot),x(id_pilot),':r',t(id_pilot),rx_pilot,'b',t(id_pilot),est_pilot,'--g');
title('TX, RX and estimated pilot');
xlabel('time /s'); ylabel('amplitude');
legend('TX pilot','RX pilot','estimated pilot');


% plot data in time domain

figure
plot(t,x,':r',t,y,'b'); hold on;
plot(t(1:length(y_clean)),y_clean,'g'); % xlim([0 8e-5]); ylim([-0.1 0.1]);
% ylim([-0.25 0.25]); xlim([0 1e-4]);
title('TX, RX and canceled data');
xlabel('time /s'); ylabel('amplitude');
legend('TX data','RX data','canceled data');



% plot_sic(x,y,y_clean,rate);
