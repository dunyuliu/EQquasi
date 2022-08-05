clear all; close all;
% This script converts results from EQquasi to rupture.dat and plots it.
% Last updated on 04/10/2021 to make the script more flexible with
% different dx and path setup. 

ht = 2; hs = 2; H = 12; 
mod = 11;
icend = 4;

[path,dx,nzz,l] = model_info(mod)

out = strcat(path,'res/catalog.dat');

header = 'header/rptheader.txt';


finalt = 0;
for ic = 1:icend
    path1 = strcat(path,'Q',num2str(ic-1),'/')
    a = load(strcat(path1,'cplot_EQquasi.txt'));
    b = load(strcat(path1,'cplot_ruptarea_trac_slip.txt'));
    % Data structure in cplot_ruptarea_trac_slip.txt;
    % Column 1   2   3   4
    %     area  slip trac_begi trac_end
    tt = load(strcat(path1,'tdyna.txt'));

    x = -l:dx:l;
    z = -l:dx:0;
    [xx,zz] =meshgrid(x,z);
    [nz,nx] = size(xx);
    ntag = 0; 
    ave_res = zeros(1,4);
    for i = 1:nx
        for j =1:nz
            area(j,i) = b((i-1)*nzz+j,1);
            slip(j,i) = b((i-1)*nzz+j,2);
            tracb(j,i) = b((i-1)*nzz+j,3);
            trace(j,i) = b((i-1)*nzz+j,4);
            if area(j,i) > 0
                ave_res(1) = ave_res(1) + area(j,i);
                ave_res(2) = ave_res(2) + slip(j,i)*area(j,i);
                ave_res(3) = ave_res(3) + tracb(j,i)*area(j,i);
                ave_res(4) = ave_res(4) + trace(j,i)*area(j,i);
            end
        end
    end
    
    ave_res(2) = ave_res(2)/ave_res(1);
    ave_res(3) = ave_res(3)/ave_res(1)/1e6;
    ave_res(4) = ave_res(4)/ave_res(1)/1e6;
    
    rec(ic,1) = ic;
    rec(ic,2) = tt(1) + finalt;
    rec(ic,3) = tt(2) + finalt;
    finalt = finalt + tt(2) + 100;
    rec(ic,4) = ave_res(1);
    rec(ic,5) = ave_res(3);
    rec(ic,6) = ave_res(4);
    rec(ic,7) = ave_res(2);
end


figure(1)
subplot(3,1,1)
plot(rec(:,1), rec(:,2), '*'); hold on;
plot(rec(:,1), rec(:,3), '*'); hold on;
title('Rupture time: start & end (s)');

subplot(3,1,2)
plot(rec(:,1), rec(:,5), '*'); hold on;
plot(rec(:,1), rec(:,6), '*'); hold on;
title('Average traction: start & end (MPa)');

subplot(3,1,3)
plot(rec(:,1), rec(:,7), '*'); hold on;
title('Average slip: start & end (m)');
xlabel ('Event Number');

delete(out);
%%output the combined file
fileID = fopen('a','w');
fprintf(fileID,'%22.14e %22.14e %22.14e %22.14e %22.14e %22.14e %22.14e\n',rec(:,:)');
fclose(fileID);
system(strcat('type'," ",'.\header\catalogheader.txt'," ",'>>'," ",out));
system(strcat('type'," ",'a'," ",'>>'," ",out));