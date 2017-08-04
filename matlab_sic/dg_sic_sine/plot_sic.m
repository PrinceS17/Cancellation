function fmax = plot_sic(x,y,y_clean,rate)

% Fourier analysis
fx = 20*log10(abs(fft(x)));
fmax = max(fx);
fx = fx - fmax;
fy = 20*log10(abs(fft(y))) - fmax;
fy_clean = 20*log10(abs(fft(y_clean))) - fmax;
% fy_clean = 20*log10(abs(fft([y_clean,zeros(1,L-length(y_clean))]))) - fmax;

% set the scale in case lengths are different
fsx = rate*(1:length(x))/length(x) - rate/2;
fsy = rate*(1:length(y))/length(y) - rate/2;
fsy_clean = rate*(1:length(y_clean))/length(y_clean) - rate/2;

% plot frequency domain
figure
plot(fsx,fftshift(fx),':r');hold on
plot(fsy,fftshift(fy),'b');xlabel('frequency /Hz');hold on
plot(fsy_clean,fftshift(fy_clean),'-g');
title('Spectrum of TX, RX and cancelled signal');
ylabel('dB');
legend('TX','RX','RX_{cancelled}');

end