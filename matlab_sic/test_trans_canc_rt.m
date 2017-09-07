% test transceiver_canceler(/_multi_tone) in real time
close all
clear all
addpath(genpath('/home/flexicon/Documents/Cancellation/Cancellation'));
build_path = {'/home/flexicon/Downloads/GNURadio/uhd/host/examples/transceiver_canceler/build';
    '/home/flexicon/Downloads/GNURadio/uhd/host/examples/transceiver_canceler_multi_tone/build';
    '/home/flexicon/Downloads/GNURadio/uhd/host/examples/tranceler_multi_tone_rt/build';
    'D:/ѧϰ/Flexicon/04code/matlab_sic/qpsk_data';
    '/mt_data'
    };
filename = {'tx_file','rx_file','y_clean_file','estimated_pilot','rx_pilot';
    'tx_file_mt','rx_file_mt','y_clean_file_mt','estimated_pilot_mt','rx_pilot_mt';
    'tx_file','rx_file','y_clean_file','estimated_pilot','rx_pilot';
    'tx_file','rx_file','y_clean_file','estimated_pilot','rx_pilot';
    'tx_file','rx_file','y_clean_file','',''
   };

id = 5;           % switch: 1 for transceiver_canceler; 2 for multi_tone; 3 for multi_tone real time; 4 for qpsk on Windows

[~,w] = size(filename);
w = 3;
a = cell(1,w);
addpath(genpath(build_path{id}));

for k = 1:w
fid = fopen(filename{id,k});              % read x, y, y_clean and pilot
a{k} = fread(fid,[1,inf],'float');
fclose(fid);
end

x0 = a{1};
y0 = a{2};
y_clean0 = a{3};
% est_pilot = a{4};
% rx_pilot = a{5};
obs = 1:6e5;

% plot pilot
rate = 2e6;   % should consistent with the C++ code!!!
f = 100e3;
spb = 20440;
L = spb;      % specify the length of a buff to visualize
st = 1;       % starting point every time

figure
h_old = cell(1,3);
num = 0;
res = 0;
ifst = 0;
while 1
    if st + L - 1 > min([length(x0),length(y0),length(y_clean0)])
        for k = 1:w
            fid = fopen(filename{id,k});              % read x, y, y_clean and pilot
            a{k} = fread(fid,[1,inf],'float');
            fclose(fid);
        end
        
        x0 = a{1};
        y0 = a{2};
        y_clean0 = a{3};
%         est_pilot = a{4};
%         rx_pilot = a{5};
        continue;
    end
    x = x0(st:st + L - 1);
    y = y0(st:st + L - 1);
    y_clean = y_clean0(st:st + L - 1);
    t = (st:st + L - 1)/rate;
    
    if id == 1
        preamble_length = 20;
        pilot_length = 400;
    else
        preamble_length = 20;
        pilot_length = 400;
    end
    
    % figure     % in fact we don't know TX's pilot here
    % id_pilot = preamble_length + 1:preamble_length + pilot_length;
    % plot(t(id_pilot),x(id_pilot),':r',t(id_pilot),rx_pilot,'b',t(id_pilot),est_pilot,'--g');
    % title('TX, RX and estimated pilot');
    % xlabel('time /s'); ylabel('amplitude');
    % legend('TX pilot','RX pilot','estimated pilot');
    
    
    % plot data in time domain
    
    % plot(t,x,':r',t,y,'b'); hold on;
    % plot(t(1:length(y_clean)),y_clean,'g'); % xlim([0 8e-5]); ylim([-0.1 0.1]);
    % % ylim([-0.25 0.25]); xlim([0 1e-4]);
    % title('TX, RX and canceled data');
    % xlabel('time /s'); ylabel('amplitude');
    % legend('TX data','RX data','canceled data');
    
    
    fc = 100e3;
    bw = 2;
    rg = 2e3;
    result = sic_db([y;y_clean],rate,fc,bw,rg),  
    if result > 30
        ifst = 1;
    end
%     if result > 55
%         break;
%     end
    if ifst
    res = res + result;
    num = num + 1;
    res/num,
    end
    
    h = plot_sic(x,y,y_clean,rate);
    for i = 1:3
        delete(h_old{i});
    end
    h_old = h;
    drawnow;
    st = st + L;
%      pause(.3);
    

end
