function y_clean = dg_sic(x,y,rate,freq,estimator_length,estimator_start,pilot_length,signal_length)
% digital cancellation for SINE wave: 
% 1, including estimate(synchronization) and cancellation
% 2, we can calculate A and A_inv outside, doen't matter in Matlab though

k = floor(estimator_length/2);
k = max(1,k);
n = pilot_length;
start = estimator_start;
L = signal_length;
N_T = round(rate/freq);
preamble_length = N_T;
preamble = x(1:N_T);

%% synchronization: detect corresponding y(n)
% Cor = xcorr(y,[zeros(1,L-n),x(start - k:start + n - k -1)]);
Cor = xcorr(y, [zeros(1, L - N_T), preamble]);
loc = pickpeaks(Cor,round(length(Cor)/2));          % ideal situation: just find peaks
loc = loc(Cor(loc) > 0.9*max(Cor));                % choose peaks that are not too small
location = min(loc);
st_rcv = location - N_T + 1;        % y(st_rcv) -> x(start)
delay = st_rcv - 1;
while st_rcv <= 0
   st_rcv = st_rcv + N_T;
end


%% estimation of pilot & figure of h
A0 = toeplitz(x(start + k - 1:start + n + k - 1),fliplr(x(start - k :start + k -1)));
A0 = fliplr(A0);
A_inv = pinv(A0);
h = A_inv*y(st_rcv:st_rcv + n)';

%% cancellation: supposed received signal - SI: A*h
N = L - max(start,st_rcv) - 2*k - n;
A = toeplitz(x(start + n + k - 1:start + N + n + k -1),fliplr(x(start + n - k:start + n + k -1)));
A = fliplr(A);
y_clean = y(st_rcv + n :st_rcv + n + N  ) - (A*h)';

%% plot the figures
figure_sic

end