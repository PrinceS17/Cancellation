clc;
clear;
close all;
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Signal Generation   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SigGen_Sync_Multitone
% SigGen creates the following variables:
% pilot_I is the in phase, pilot portion of the transmitted signal
% pilot_Q is the quadrature, pilot portion of the transmitted signal
% sig_I is the in phase transmitted signal
% sig_Q is the quadrature transmitted signal
% N_tot_pilot is the length of the pilot signal
% N_tot_signal si the length of the actual signal
% N-tot = N_tot_pilot + N_tot_signal
% vout is the in phase waveform on the receiver side
%% don't know what the vout exactly is

IQ_Rate = 2E6; % Sampling rate on both the tx side and rx side


f1 = 400E3; % tone frequency 1
f2 = 500E3; % tone frequency 2

T = 50000; % period for the two tone signal
T_N_tot = N_tot/IQ_Rate; % total time duration for pilot and signal
T_N_tot_pilot = (N_tot/IQ_Rate)*N_tot_pilot/N_tot; % time duration for pilot
%% T_N_tot_pilot = N_tot_pilot/IQ_Rate;?

T_shift = 0; % softwared-defined delay between the triggering instant and the first data points (ns)
% T_shift needs to be smaller than 2500ns=2.5us
Sample_shift = T_shift/0.25; % number of samples that the original sequence 'vout' shifted
L_vout = length(vout);
for i=1:(L_vout-Sample_shift)
    vout_shift(i) = vout(i+Sample_shift);
end
%% vout_shift = vout(Sample_shift + 1:end);

 
% Find the averaged pulse response from the single frame in 2M Sa/s rate
% Length of one period for combined 400 kHz and 500 kHz tones =
% 1/gcd(400k,500k) = 1E-5
% # of samples per period = 2M Sa/s * 1E-5s = 20
% N_period = IQ_Rate*T; % Number of samples per period 
% Num_Cyc = T_N_tot_pilot/(T); % number of cycles in the pilot tone 0.00025s/1E-5s = 25 
% AVG = 0.8*Num_Cyc; % 80% of the total number of cycles
% v_avg_temp = 0;
% v_avg = zeros(1,N_period);
% 
%  for p=1:N_tot_pilot
%      for q=1:(AVG-1) % the number of average is AVG-1, as the entire sequence is advanced by Sample_shift
%          v_avg_temp = vout_shift((q+Num_Cyc-AVG)*N_period+p)+ v_avg_temp;
%      end
%      v_avg(p) = v_avg_temp/AVG;
%      v_avg_temp = 0;
%  end
 
% name v_pilot_avg as the averaged pilot tone at 2M Sa/s rate
v_pilot_avg = vout(1:N_tot_pilot);

% Separate the two-tone signal from the pilot in a frame capture
v_sig = vout_shift(N_tot_pilot+1:N_tot);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Decimate the captured signal from scope to the desired rate %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMOVED: This section is unnecessary because the IQ Rate for Rx and Tx are both
% 2MSa/s

% DwSaRt = 100; % Down sampling rate from the original 4G Sa/s rate
% 
% % Down-sampled pilot signal v_pilot_avg_d
% for i=1:(length(v_pilot_avg)/DwSaRt)
%     v_pilot_avg_d(i) = v_pilot_avg(1+(i-1)*DwSaRt);
% end
% 
% % Down-sampled two-tone RX signal v_sig_d
% for i=1:(length(v_sig)/DwSaRt)
%     v_sig_d(i) = v_sig(1+(i-1)*DwSaRt);
% end

v_pilot_avg_d = v_pilot_avg;
v_sig_d = v_sig; 

figure(1)
stem(v_pilot_avg_d)
xlabel('Samples (-)')
ylabel('Sample Magnitude (-)')
title('Pilot tone sampled at 5 MSa/s');


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Nonlinear Coefficient   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
Delay_L_est = 42; % Estimated length of the delay
Nonlinear_dimen = 4; % number of nonlinear terms that taken into account.
                     % 4 means 1st, 2nd, 3rd, and 4th order nonlinearity are
                     % considered.
ofst = 50; % offset when look at the input sequence. 
            % ofst needs to be larger than Delay_L_est
% Jin's
% find the matrix that acts on the nonlinear coefficient  
X = zeros(length(pilot_I)-ofst,(Nonlinear_dimen-1)*(Delay_L_est+1)+Delay_L_est);
for m=1:(length(pilot_I)-ofst)
    for p=1:Nonlinear_dimen
        for n=1:(Delay_L_est+1)
            X(m,(p-1)*(Delay_L_est+1)+n) = (pilot_I(ofst+m-(n-1)))^(1+(p-1));
        end
    end
end

%% in the paper only odd order terms are considered? here is different?

%% Tingjun's
% Delay_L_est = 150;
% Nonlinear_dimen_vals = [1,2,3,5,6,7,8,9,10,11];
% Nonlinear_dimen = length(Nonlinear_dimen_vals);
% Nonlinear_dimen_index = [1:1:Nonlinear_dimen];
% ofst = 200;
% 
% X = zeros(length(pilot_I)-ofst,(Nonlinear_dimen-1)*(Delay_L_est+1)+Delay_L_est);
% % for m=1:(length(pilot_I)-ofst)
% for m=1:Delay_L_est
%     for p=Nonlinear_dimen_index
%         for n=1:(Delay_L_est+1)
%             X(m,(p-1)*(Delay_L_est+1)+n) = (pilot_I(ofst+m-(n-1))).^(1+(Nonlinear_dimen_vals(p)-1));
%         end
%     end
% end

% end of Tingjun

% find the nonlinear coefficient using pseduo-inverse of matrix X
pseduo_inv_X = pinv(X);
h = zeros(1,(Delay_L_est+1)*Nonlinear_dimen);
% find the output sequence that corresponds to the above matrix
for m=1:((Delay_L_est+1)*Nonlinear_dimen)
    for n=1:(length(pilot_I)-ofst)
        h(m) = v_pilot_avg_d(ofst+n)*pseduo_inv_X(m,n)+h(m);
    end
end

h = h';

%% why? just h = pseduo_inv_X(1:....,:)*v_pilot_avg_d(ofst+1:length(pilot_I));


% find the output sequence that corresponds to the above matrix
L_v_pilot_avg_d = length(v_pilot_avg_d);
v_pilot_avg_d_shifted = v_pilot_avg_d(ofst+1:L_v_pilot_avg_d);

% find the output signal based on the nonlinear coefficient
v_pilot_calc = X*h;

 figure(2)
 stem(v_pilot_avg_d_shifted)
 hold on;
 stem(v_pilot_calc)
 xlabel('Samples (-)')
 ylabel('Sample Magnitude (-)')
 legend('meas. pilot signal','calc. pilot using nonlinear coefficient')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Perform Digital SIC    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_tot = 500; % length of the sequence that is taken into account
ofst_sig = 250;
tic

A_TT = zeros(N_tot,(Nonlinear_dimen-1)*(Delay_L_est+1)+Delay_L_est+1);
for m=1:N_tot
    for p=1:Nonlinear_dimen
        for n=1:(Delay_L_est+1)
            A_TT(m,(p-1)*(Delay_L_est+1)+n) = (sig_I(ofst_sig+m-(n-1)))^(1+(p-1));
        end
    end
end


v_sig_d_hat = A_TT*h; % cancellation signal
toc

% find the output sequence that corresponds to the above matrix
v_sig_d_shifted = zeros(1,N_tot);
for m=1:N_tot
    v_sig_d_shifted(m) = v_sig_d(ofst_sig+m);
end
v_sig_d_shifted = v_sig_d_shifted';

figure(3)
stem(v_sig_d_shifted)
hold on;
stem(v_sig_d_hat,'r')
xlabel('Samples (-)')
ylabel('Sample Magnitude (-)')
title('SI and Cancellation Signal, 40M Sa/second');
legend('SI','Cancellation Signal')

% find the original sequence that corresponses to the above matrix
sig_I_shifted = sig_I(ofst_sig+1:ofst_sig+N_tot);
sig_I_shifted = sig_I_shifted';

% residual SI
Re_SI = v_sig_d_shifted-v_sig_d_hat;
Re_SI = Re_SI';

% FFT of the output signal spectrum
f = 0:1:(N_tot-1);
f = IQ_Rate/N_tot.*f;
% Original Two-tone
V_sig_I_shifted = fft(sig_I_shifted)/length(f)*2;
mag_V_sig_I_shifted= abs(V_sig_I_shifted);
% SI
V_v_sigd_shifted = fft(v_sig_d_shifted)/length(f)*2;
mag_V_v_sigd_shifted = abs(V_v_sigd_shifted);
% Residual SI
V_Re_SI = fft(Re_SI)/length(f)*2;
mag_V_Re_SI = abs(V_Re_SI);
% Calculate the power from the voltage magnitude
P_sig_I_shifted = 10*log10(mag_V_sig_I_shifted.^2/100*1000);
P_v_sigd_shifted = 10*log10(mag_V_v_sigd_shifted.^2/100*1000);
P_v_sigd_shifted = P_v_sigd_shifted';
P_Re_SI = 10*log10(mag_V_Re_SI.^2/100*1000);
% plot the power versus frequency

figure(4)
plot(f,P_v_sigd_shifted)
hold on;
grid on;
plot(f,P_Re_SI,'r--','linewidth',2)
xlim([0 f2*4])
xlabel('Frequency (Hz)')
ylabel('Power (dBm)')
title('Received Signal with and without digital SIC');
legend('before SIC','residual SI after digital SIC')

% figure(5)
% stem(h)
% hold on;
% grid on;
% xlabel('Number of nonlinear coefficent (-)')
% ylabel('Magnitude (-)')
% title('Volterra Series Nonlinear Coefficient');

% Calculate the SIC at the main tone and IM3 tone frequencies
SIC_dB = P_v_sigd_shifted - P_Re_SI;
index_SI_main1 = find(f==1.6e6);
index_SI_main2 = find(f==2e6);
index_SI_IM3_1 = find(f==1.2e6);
index_SI_IM3_2 = find(f==2.4e6)
SIC_SI_main1 = SIC_dB(index_SI_main1)
SIC_SI_main2 = SIC_dB(index_SI_main2)
SIC_SI_IM3_1 = SIC_dB(index_SI_IM3_1)
SIC_SI_IM3_2 = SIC_dB(index_SI_IM3_2)

% output peak
max(P_Re_SI)

