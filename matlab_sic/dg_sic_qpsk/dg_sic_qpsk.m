function [y_clean, MSE] = dg_sic_qpsk(x, y, rate, N_T, preamble, pilot, data_length, estimator_length, start)
% digital cancellation for QPSK, start is where we begin to estimate the
% channel, x, y should be a later longer than the buffer

preamble_length = length(preamble);
k = floor(estimator_length/2);
lx = length(x);
ly = length(y);
n = length(pilot) + length(preamble) - start - k;    % use pilot after start to estimate h 

% synchronization: actually shouldn't deviate right location farther than
% length of estimator, or we can't get the right channel!!

Cor = xcorr((y), [zeros(1,ly - length(preamble)),preamble]);      % for complex number
% location = pickpeaks(abs(Cor),1,0);                             % find the max's index, hard judge?

loc = pickpeaks(Cor,round(length(Cor)/4));          % ideal situation: just find peaks
loc = loc(Cor(loc) > 0.9*max(Cor));                % choose peaks that are not too small
location = min(loc);

delay = location - length(preamble);          % calculate the delay
st_rcv = start + delay;                       % get the start of pilot in received signal
delay,
st_rcv,
% st_rcv = start ;

% % plot the preamble
% figure
% t = (start:length(preamble))/samp_rate;
% plot(t, preamble(start:end)); hold on;
% plot(t, y(st_rcv:st_rcv + length(preamble) - start), ':r');


% channel estimation
A0 = toeplitz(x(start + k - 1:start + n + k - 1),fliplr(x(start - k :start + k -1)));
A0 = fliplr(A0);
A_inv = pinv(A0);
h = A_inv*y(st_rcv:st_rcv + n)';
% self_test = y(st_rcv:st_rcv + n) - conj(A*h)';   % ' is conjugate transpose !!


% digital cancellation
% N = ly - max(st_rcv,start) - 2*k - n;

N = data_length;        % do the cancellation from st_rcv:st_rcv + data_length
A = toeplitz(x(start + n + k - 1:start + N + n + k -1),fliplr(x(start + n - k:start + n + k -1)));
A = fliplr(A);
y_clean = y(st_rcv + n:st_rcv + n + N) - (A*h)';
MSE = norm(y_clean)^2/(N + 1);

% plot the figures
figure_sic

end