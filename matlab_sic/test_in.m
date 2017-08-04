fid = fopen('qpsk_output');              % read y_clean
a = fread(fid,[2,inf],'float');
fclose(fid);
plot(a(1,:));