function result = sic_db(signal, rate, fc, bw)
% calculate cancellation result: signal is matrix of RX and y_clean, rate is % sampling rate, 
% fc is center frequency, bw is bandwidth which can be calculated theoretically, 
% result is column vector of cancellation results

P_db = zeros(2,1);
for i = 1:2
    x = signal(i,:);
    fx = abs(fft(x))/length(x)*2;
    Px = 10*log10(fx.^2/100*1000); 
    Px = fftshift(Px);
    
    fl_id = ((fc - bw/2)/rate + 0.5)*length(x);
    fr_id = ((fc + bw/2)/rate + 0.5)*length(x);
    
    P_db(i) = mean(Px(fl_id:fr_id));
end
result = P_db(1) - P_db(2);

end
    
