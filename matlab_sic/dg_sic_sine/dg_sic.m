function y_clean = dg_sic(x, y, rate, freq, data_length, estimator_length, estimator_start, pilot_length, signal_length, dim)
% digital cancellation for SINE wave: 
% 1, including estimate(synchronization) and cancellation
% 2, we can calculate A and A_inv outside, doen't matter in Matlab though
% 3, now include nonlinear cancellation

if nargin < 10
    dim = 1;            % use linear cancellation by default
end

k = floor(estimator_length/2);
k = max(1,k);
lx = length(x);
ly = length(y);
n = pilot_length;
start = estimator_start;
N_T = round(rate/freq);
preamble_length = N_T;
preamble = x(1:N_T);

%% synchronization: detect corresponding y(n)
% Cor = xcorr(y,[zeros(1,L-n),x(start - k:start + n - k -1)]);
Cor = xcorr(y, [zeros(1, ly - N_T), preamble]);
loc = pickpeaks(Cor,round(length(Cor)/2));         % ideal situation: just find peaks
loc = loc(Cor(loc) > 0.9*max(Cor));                % choose peaks that are not too small
location = min(loc);
delay = location ;                           % y(st_rcv) -> x(start)
if delay < 0
    'Error: delay is negative!'
end

% delay = k - 1;
st_rcv = start + delay;

%% estimation of pilot & figure of h
dim,
A0 = mat_generation(x, start, n, k, dim);  
A_inv = pinv(A0);
h = A_inv*y(st_rcv:st_rcv + n)';

%% output for debug
rank_A = rank(A0),
norm_h = norm(h),

%% cancellation: supposed received signal - SI: A*h
N = data_length;
A = mat_generation(x, start + n, N, k, dim);
y_clean = y(st_rcv + n :st_rcv + n + N  ) - (A*h)';

%% plot the figures
figure_sic

end