% plot TX & RX figures
% needed variables: h, rate, start, n (pilot or pilot + preamble), A0,
% N(data length), x, y, A, y_clean

%% plot h
figure
stem(h); title('coefficient of h'); 

%% plot pilot and estimated pilot
t = (start:start + n)/rate;
figure
plot(t,x(start:start + n),':r',t,y(st_rcv:st_rcv + n),'b',t,(A0*h)','--g');
title('TX & RX pilot and estimated pilot');
xlabel('time /s'); ylabel('amplitude');
legend('TX pilot','RX pilot','estimated pilot');

%% plot TX & RX data
t1 = (start + n:start + n + N)/rate;        % N is the cancelled data length 
figure 
plot(t1,x(start + n:start + n + N),':r',t1,y(st_rcv + n:st_rcv + n + N),'b',t1,(A*h)','--g');
title('TX, RX and estimated data');
xlabel('time /s'); ylabel('amplitude');
legend('TX data','RX data','estimated data');

%% plot cancellation result: frequency domain
plot_sic(x,y,y_clean,rate);