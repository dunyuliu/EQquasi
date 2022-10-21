clear all; close all;

x = -35:1:35;
z = -25:1:0;

[xx,zz] = meshgrid(x,z);

[n,m] = size(xx);

for i = 1:n
    for j = 1:m
        xcoor = xx(i,j);
        zcoor = zz(i,j);
        
        b(i,j) = 0.03;
        
        if abs(zcoor)>=18 || abs(xcoor)>=32 || abs(zcoor)<=2.0
            a(i,j) = 0.04;
        elseif (abs(zcoor)<=16 && abs(zcoor)>=4 && abs(xcoor)<=30)
            a(i,j) = 0.004;
        else
            tmp1 = abs(abs(zcoor)-2-2-12/2) - 12/2;
            tmp1 = tmp1/2;
            tmp2 = (abs(xcoor)-30)/2;
            a(i,j) = 0.004 + max(tmp1,tmp2)*(0.04-0.004);
        end
    end
end

figure(1)
contourf(xx,zz,a-b); colorbar; axis equal;