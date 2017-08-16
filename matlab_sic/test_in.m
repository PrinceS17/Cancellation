path = '/home/flexicon/Downloads/GNURadio/uhd/host/examples/qpsk_tranceler/build/';
fid = fopen(strcat(path,'tx_file'));              % read tx
a1 = fread(fid,[1,inf],'float');
fclose(fid);

fid = fopen(strcat(path,'rx_file'));              % read rx
a2 = fread(fid,[1,inf],'float');
fclose(fid); 

t = 1:length(a1);
figure
plot(t,a1,'--r',t,a2,'b');