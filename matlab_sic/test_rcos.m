% test with raised cosine filter
%% 1, use rcosdesign
symbol = [1 + 1i, 1 - 1i, -1 + 1i, -1 - 1i];
ls = 1e3;
sps = 16;
beta = 0.5;
span = 8;

x = symbol(randi(4,1,ls));
rcf = rcosdesign(beta,span,sps,'normal');
lf = length(rcf);

% tx = upsconv1(x,rcf,1e3);
x1 = upsample(x,sps);
tx = conv(x1,rcf);
tx = tx(1:length(tx) - (lf - 1));    % delete lf - 1 pts at all

%% 2, use qpsk_generation
tx2 = qpsk_generation(x,beta,span,sps);

% plot 
figure
L = sps*ls;
subplot(2,1,1),plot(1:L,real(tx),'--r',1:L,real(tx2),':b');
subplot(2,1,2),plot(-L/2 + 1:L/2,fftshift(abs(fft(tx))),'--r',-L/2 + 1:L/2,fftshift(abs(fft(tx2))),':b');
