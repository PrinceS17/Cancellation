function y_clean = dg_sic(x,y,rate,freq,estimator_length,estimator_start,pilot_length,signal_length)
% digital cancellation for SINE wave: 
% 1, including estimate(synchronization) and cancellation
% 2, we can calculate A and A_inv outside, doen't matter in Matlab though

k = floor(estimator_length/2);
n = pilot_length;
start = estimator_start;
L = signal_length;
N_T = round(rate/freq);
preamble = x(start:start + N_T - 1);


%% synchronization: detect corresponding y(n)
% Cor = xcorr(y,[zeros(1,L-n),x(start - k:start + n - k -1)]);
Cor = xcorr(y, [zeros(1, L - N_T), preamble]);
loc = pickpeaks(Cor,round(length(Cor)/2));          % ideal situation: just find peaks
loc = loc(Cor(loc) > 0.85*max(Cor));                % choose peaks that are not too small
location = min(loc);
st_rcv = location - n + k + 1;        % y(st_rcv) -> x(start)
if st_rcv <= 0
    st_rcv = st_rcv + n;
end

%% estimation of pilot & figure of h
A = toeplitz(x(start + k - 1:start + n + k - 1),fliplr(x(start - k :start + k -1)));
A = fliplr(A);
A_inv = pinv(A);
h = A_inv*y(st_rcv:st_rcv + n)';
figure
stem(h); title('coefficient of h'); 

%% plot pilot and estimated pilot
t = (start:start + n)/rate;
figure
plot(t,x(start:start + n),':r',t,y(st_rcv:st_rcv + n),'b',t,(A*h)','--g');
title('TX & RX pilot and estimated pilot');
xlabel('time /s'); ylabel('amplitude');
legend('TX pilot','RX pilot','estimated pilot');

%% cancellation: supposed received signal - SI: A*h
N = L - max(start,st_rcv) - 2*k - n;
A = toeplitz(x(start + n + k - 1:start + N + n + k -1),fliplr(x(start + n - k:start + n + k -1)));
A = fliplr(A);
y_clean = y(st_rcv + n :st_rcv + n + N  ) - (A*h)';

%% plot TX & RX data
t1 = (start + n:start + n + N)/rate;
figure 
plot(t1,x(start + n:start + n + N),':r',t1,y(st_rcv + n:st_rcv + n + N),'b',t1,(A*h)','--g');
title('TX, RX and estimated data');
xlabel('time /s'); ylabel('amplitude');
legend('TX data','RX data','estimated data');


end