fid = fopen('sine_wave_2M_g');            % read y
a1 = fread(fid,[1,inf],'float');
a1 = a1(2000000:2020000);
fclose(fid);
fid = fopen('sine_wave_2M_new','w');
fwrite(fid,a1,'float');
fclose(fid);