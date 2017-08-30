function result = sic_db(signal, rate, fc, bw, rg)
% calculate cancellation result: signal is matrix of RX and y_clean, rate is % sampling rate, 
% fc is center frequency, bw is bandwidth which can be calculated theoretically, 
% result is column vector of cancellation results; rg is the search range
% for true fc (may be a shift off given fc)

P_db = zeros(2,1);
[~, L] = size(signal);
fl_id = round(((fc - rg/2)/rate + 0.5)*L);
fr_id = round(((fc + rg/2)/rate + 0.5)*L);
for i = 1:2
    x = signal(i,:);
    fx = abs(fft(x))/L*2;
    Px = 10*log10(fx.^2/100*1000); 
    Px = fftshift(Px);
    
     % search from [-bw,0] to [0,bw]
    
    if i == 1
        [~,peak_id] = max(Px(fl_id:fr_id));             % for searching max
        k = round(bw * L /2 /rate);
        fr_id = fl_id + peak_id - 1 + k;
        fl_id = fl_id + peak_id - 1 - k;
       
    end
    
    P_db(i) = mean(Px(fl_id:fr_id));
end
result = P_db(1) - P_db(2);

end
    
