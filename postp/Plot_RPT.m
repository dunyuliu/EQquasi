clear all; close all;
% This script converts results from EQquasi to rupture.dat and plots it.
% Last updated on 04/10/2021 to make the script more flexible with
% different dx and path setup. 

ht = 2; hs = 2; H = 12; 
mod = 31;
ic = 1;

[path,dx,nzz,l] = model_info(mod)

out = strcat(path,'res/rupture.dat');
header = 'header/rptheader.txt';

path1 = strcat(path,'Q',num2str(ic-1),'/')
a = load(strcat(path1,'cplot_EQquasi.txt'));
tt = load(strcat(path1,'tdyna.txt'));

x = -l:dx:l;
z = -l:dx:0;
zcenter = 0;
% if mod == 12
%     z=-2*l:dx:0;
%     zcenter = -60;
% end
[xx,zz] =meshgrid(x,z);
[nz,nx] = size(xx);
ntag = 0; 
for i = 1:nx
    for j =1:nz
        rpt(j,i) = a((i-1)*nzz+j,16)-tt(1,1);
        if rpt(j,i)<0
            rpt(j,i) = 100000;
        end
        if (abs(xx(j,i))<= l/2+2*ht && abs(zz(j,i)-zcenter)<= H +2*ht +hs )
            ntag = ntag + 1;
            rptout(ntag,1) = xx(j,i)*1e3;
            rptout(ntag,2) = -(zz(j,i))*1e3;
            rptout(ntag,3) = rpt(j,i);
        end
    end
end

ll = 0:10:300;
figure(1)
contour(xx,zz,rpt,ll, 'ShowText' ,'on');colorbar;axis equal;hold on;

delete(out);
%%output the combined file
fileID = fopen('a','w');
fprintf(fileID,'%22.14e %22.14e %22.14e \n',rptout(:,:)');
fclose(fileID);
system(strcat('type'," ",'.\header\rptheader.txt'," ",'>>'," ",out));
system(strcat('type'," ",'a'," ",'>>'," ",out));