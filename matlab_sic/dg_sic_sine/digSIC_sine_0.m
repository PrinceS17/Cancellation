% digital cancellation for SINE wave

%% constant & signal
symbol_num = 50;                     % num of symbols
f = 400e3;                           % wave frequency
rate = 10e6;                         % sampling rate
T = 1/f;                             % period of a symbol
phi = pi*0.6;                        % phase of sine wave
SNR = 20;
L = round(symbol_num*T*rate);        % length of the whole signal
t = (1:L)/rate;
x = 0.3*cos(2*pi*t/T + phi);         % TX signal
y = 0.1*x + 0.005*[0,x(1:end-1)] - 0.002*[x(2:end),0];    % generated for simulation
y = awgn(y,SNR,'measured');                               % add noise (measured is right)
delay = 1e2;                                              % delay between x and y
y = [y(delay + 1:end),0.001*randn(1,delay)];              % receiver delay 


%% estimation: calculate response h
k = 2;                              % number of pts before and after current, 2k is estimator length 
start = 1.5*k;                      % when the pilot begin 
n = round(10*rate*T);               % pilot length
A = toeplitz(x(start + k - 1:start + n + k - 1),fliplr(x(start - k :start + k -1)));
A = fliplr(A);
A_inv = pinv(A);

Cor = xcorr(y,[zeros(1,L-n),x(start - k:start + n - k -1)]);    % detect corresponding y(n)
loc = pickpeaks(Cor,round(length(Cor)/2));          % ideal situation: just find peaks
loc = loc(Cor(loc) > 0.85*max(Cor));                % choose peaks that are not too small
location = min(loc);
st_rcv = location - n + k + 1;                                    % y(st_rcv) matches x(start)
if st_rcv <= 0
    st_rcv = st_rcv + n;
end
h = A_inv*y(st_rcv:st_rcv + n)';


%% cancellation: supposed received signal - SI: A*h
lend = L/10;
N = L - max(start,st_rcv) - 2*k - n - lend;                         % the irregular end will influence the cancellation
A = toeplitz(x(start + n + k - 1:start + N + n + k -1),fliplr(x(start + n - k:start + n + k -1)));
A = fliplr(A);
y_clean = y(st_rcv + n:st_rcv + n + N  ) - (A*h)';      % do the cancellation after pilot


%% plot: time domain
subplot(2,1,1),plot(t,x,':r');xlabel('time /s');hold on
plot(t,y,'b');hold on
plot(t(start + n:start + n + N),y_clean,'-g');
legend('input x','output y','y_{clean}');


%% plot: Fourier analysis
fx = 20*log10(abs(fft(x)));             % use dB unit
fmax = max(fx);
fx = fx - fmax;                         % normalize the scale
fy = 20*log10(abs(fft(y))) - fmax;
fy_clean = 20*log10(abs(fft([y_clean,zeros(1,L-length(y_clean))]))) - fmax;
fscale = rate*(1:L)/L - rate/2;
subplot(2,1,2),plot(fscale,fftshift(fx),':r');hold on
plot(fscale,fftshift(fy),'b');xlabel('frequency /Hz');hold on
plot(fscale,fftshift(fy_clean),'-g');
ylim([-100,10]);ylabel('dB');
legend('FFT\{x\}','FFT\{y\}','FFT\{y_{clean}\}');

% figure(2);
% f1 =  rate*(1:2*k)/(2*k) - rate/2;
% plot(f1,fftshift(abs(fft(h))));

