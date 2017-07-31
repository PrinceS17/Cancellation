function y_clean = dg_sic(x,y,estimator_length,estimator_start,pilot_length,signal_length)
% digital cancellation for SINE wave: 
% 1, including estimate(synchronization) and cancellation
% 2, we can calculate A and A_inv outside, doen't matter in Matlab though

k = floor(estimator_length/2);
n = pilot_length;
start = estimator_start;
L = signal_length;

A = toeplitz(x(start + k - 1:start + n + k - 1),fliplr(x(start - k :start + k -1)));
A = fliplr(A);
A_inv = pinv(A);

% synchronization: detect corresponding y(n)
Cor = xcorr(y,[zeros(1,L-n),x(start - k:start + n - k -1)]);
loc = pickpeaks(Cor,round(length(Cor)/2));          % ideal situation: just find peaks
loc = loc(Cor(loc) > 0.85*max(Cor));                % choose peaks that are not too small
location = min(loc);
st_rcv = location - n + k + 1;        % y(st_rcv) -> x(start)
if st_rcv <= 0
    st_rcv = st_rcv + n;
end
h = A_inv*y(st_rcv:st_rcv + n)';

% cancellation: supposed received signal - SI: A*h
N = L - max(start,st_rcv) - 2*k - n;
A = toeplitz(x(start + n + k - 1:start + N + n + k -1),fliplr(x(start + n - k:start + n + k -1)));
A = fliplr(A);
y_clean = y(st_rcv + n :st_rcv + n + N  ) - (A*h)';


end
