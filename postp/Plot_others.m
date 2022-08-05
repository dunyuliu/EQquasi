clear all; close all;
% This script converts results from EQquasi to rupture.dat and plots it.
% Last updated on 04/10/2021 to make the script more flexible with
% different dx and path setup. 

l = 60; ht = 2; hs = 2; H = 12; 
mod = 4;
ic = 1;
col = 3;
[path,dx,nzz,l] = model_info(mod)
out = strcat(path,'rupture.dat');


path1 = strcat(path,'Q',num2str(ic-1),'/')
a = load(strcat(path1,'cplot_EQquasi.txt'));
tt = load(strcat(path1,'tdyna.txt'));

x = -l:dx:l;
z = -l:dx:0;
[xx,zz] =meshgrid(x,z);
[nz,nx] = size(xx);
ntag = 0; 
for i = 1:nx
    for j =1:nz
        rpt(j,i) = a((i-1)*nzz+j,col);%-tt(1,1);
        %if rpt(j,i)<0
        %    rpt(j,i) =1d9;
        %end
        if (abs(xx(j,i))<= l/2+2*ht && abs(zz(j,i))<= H +2*ht +hs )
            ntag = ntag + 1;
            rptout(ntag,1) = xx(j,i)*1e3;
            rptout(ntag,2) = -(zz(j,i))*1e3;
            rptout(ntag,3) = rpt(j,i);
        end
    end
end

ll = 0:2:300;
figure(1)
%contour(xx,zz,rpt,ll);colorbar;axis equal; 
contourf(xx,zz,rpt);colorbar;axis equal;

delete rupture.dat
%%output the combined file
fileID = fopen('a','w');
fprintf(fileID,'%22.14e %22.14e %22.14e \n',rptout(:,:)');
fclose(fileID);
system(strcat('type'," ",'rptheader.txt'," ",'>>'," ",out));
system(strcat('type'," ",'a'," ",'>>'," ",out));