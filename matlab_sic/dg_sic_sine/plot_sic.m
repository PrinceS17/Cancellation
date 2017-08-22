function h = plot_sic(x,y,y_clean,rate)

% Fourier analysis
fx = abs(fft(x))/length(x)*2;
fy = abs(fft(y))/length(y)*2;
fy_clean = abs(fft(y_clean))/length(y_clean)*2;     % normalize: sum of whole spectrum is 1; *2: one side spectrum
    
Px = 10*log10(fx.^2/100*1000);                      % P = U^2/2R *1000 (mW)
Py = 10*log10(fy.^2/100*1000);
Py_clean = 10*log10(fy_clean.^2/100*1000);

% set the scale in case lengths are different
fsx = rate*(1:length(x))/length(x) - rate/2;
fsy = rate*(1:length(y))/length(y) - rate/2;
fsy_clean = rate*(1:length(y_clean))/length(y_clean) - rate/2;

% plot frequency domain
% figure    % comment only when using test_trans_canc_rt
h1 = plot(fsx,fftshift(Px),':r');hold on
h2 = plot(fsy,fftshift(Py),'b');xlabel('frequency /Hz');hold on
h3 = plot(fsy_clean,fftshift(Py_clean),'-g');
title('Spectrum of TX, RX and cancelled signal');

ylim([-120 0]);
ylabel('dB');
legend('TX','RX','RX_{cancelled}');
h = {h1,h2,h3};

end