clear all; close all;
out = 'rpt.txt';
ic = 2;
path = strcat('../work-ksi0.025-dt0.05/C',num2str(ic-1),'/')
a = load(strcat(path,'cplot_EQquasi.txt'));
tt = load(strcat(path,'tdyna.txt'));
dx = 1;
x = -60:dx:60;
z = -60:dx:0;
[xx,zz] =meshgrid(x,z);
[nz,nx] = size(xx);
ntag = 0; 
for i = 1:nx
    for j =1:nz
        rpt(j,i) = a((i-1)*61+j,16)-tt(1,1);
        if rpt(j,i)<0
            rpt(j,i) =300;
        end
       
        ntag = ntag + 1;
        rptout(ntag,1) = xx(j,i)*1e3;
        rptout(ntag,2) = -(zz(j,i))*1e3;
        rptout(ntag,3) = rpt(j,i);

    end
end

ll = 0:0.5:300;
figure(1)
contour(xx,zz,rpt,ll);colorbar;axis equal; 

%%output the combined file
system(strcat('rm'," ",out));
fileID = fopen('tmp','w');
fprintf(fileID,'%22.14e %15.7e %15.7e\n',rptout(:,1:3));
fclose(fileID);
system(strcat('touch'," ",out));
system(strcat('cat'," ",'rptheader.txt'," ",'>>'," ",out));
system(strcat('cat'," ",'tmp'," ",'>>'," ",out));