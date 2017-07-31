% test tranceiver.cpp
addpath(genpath('C:\Users\sjh\Documents\Visual Studio 2010\Projects\dg_sic_1.3\dg_sic_1.3'));
fid = fopen('txout');               % read x
a0 = fread(fid,[1,inf],'float');
fclose(fid);

fid = fopen('gen_data1');            % read y
a1 = fread(fid,[1,inf],'float');
fclose(fid);

fid = fopen('output');              % read y_clean
a2 = fread(fid,[1,inf],'float');
fclose(fid);
% subplot(2,1,1);
% plot(a1(1,:),':r'); hold on
% plot(a2(1,:));
% legend('y','y_{clean}');

% plot
rate = 2e6;   % should consistent with the C++ code!!!
f = 100e3;
L = length(a1);
t = (1:L)/rate;
% x = 0.5*sin(2*pi*f*t);
x = a0;
y = a1;
y_clean = a2;
plot_sic(x,y,y_clean,t,1,L,rate);
