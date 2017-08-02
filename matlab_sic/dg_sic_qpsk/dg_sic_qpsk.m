function [y_clean, MSE] = dg_sic_qpsk(x, y, preamble, pilot, estimator_length, start)
% digital cancellation for QPSK, start is where we begin to estimate the
% channel

k = floor(estimator_length/2);
lx = length(x);
ly = length(y);
n = length(pilot) + length(preamble) - start - k;    % use pilot after start to estimate h 

% synchronization: actually shouldn't deviate right location farther than
% length of estimator, or we can't get the right channel!!

Cor = xcorr(conj(y), [zeros(1,ly - length(preamble)),preamble]);      % for complex number
% location = pickpeaks(abs(Cor),1,0);                                   % find the max's index, hard judge?

loc = pickpeaks(Cor,round(length(Cor)/4));          % ideal situation: just find peaks
loc = loc(Cor(loc) > 0.9*max(Cor));                % choose peaks that are not too small
location = min(loc);

delay = location - length(preamble);          % calculate the delay
st_rcv = start + delay;                       % get the start of pilot in received signal
% st_rcv = start ;

% channel estimation
A = toeplitz(x(start + k - 1:start + n + k - 1),fliplr(x(start - k :start + k -1)));
A = fliplr(A);
A_inv = pinv(A);
h = A_inv*conj(y(st_rcv:st_rcv + n)');
% self_test = y(st_rcv:st_rcv + n) - conj(A*h)';   % ' is conjugate transpose !!


% digital cancellation
N = ly - st_rcv - 2*k - n;
A = toeplitz(x(start + n + k - 1:start + N + n + k -1),fliplr(x(start + n - k:start + n + k -1)));
A = fliplr(A);
y_clean = y(st_rcv + n:st_rcv + n + N) - conj(A*h)';
MSE = norm(y_clean)^2/(N + 1);

end