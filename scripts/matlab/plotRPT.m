function plotRPT(bp, mod)
    txt = 'Generating rupture.dat ... ...'
    % This script converts results from EQquasi to rupture.dat and plots it.
    % Change log:
    % 20220812: add plots for other components.
    % 20210410: make the script more flexible with different dx and path setup. 

    % Data structure of cplot_EQquasi.txt.
    % Col:
    % 1, 2, 3,   4,     5, 6, 7, 8, 9, 
    % x, z, slr, theta, ts,td,tn

    %bp   = 7; % benchmark problem id. 
    %mod  = 3;
    ic   = 1;
    col  = 16; % Use the col_th column data of cplot_EQquasi.txt. 

    if     bp == 5 
        [path,dx,nzz,l] = model_info_bp5(mod);
        ht   = 2; hs = 2; H = 12; 
        x    = -l:dx:l;
        z    = -l:dx:0;
        zcenter = 0;
        L    = l/2;
        ll = 0:10:300; % contour numbers.
    elseif bp == 7
        [path,dx,nzz,l] = model_info_bp7(mod);
        ht   = 0; hs = 0; 
        x    = -l/2:dx:l/2;
        z    = -l/2:dx:l/2;
        zcenter = 0;   
        L    = 400/1000;
        H    = L;
        ll   = 0:0.1:3;
    end

    out      = strcat(path,'rupture.dat');

    path1    = strcat(path,'Q',num2str(ic-1),'/');
    a        = load(strcat(path1,'cplot_EQquasi.txt'));
    tt       = load(strcat(path1,'tdyna.txt'));

    [xx,zz]  = meshgrid(x,z);
    [nz,nx]  = size(xx);

    ntag     = 0; 
    ntag1    = 0; % for bp7. Because the output of EQquasi has a larger region than required.
    for i = 1:nx
        for j =1:nz
            if col == 16
                rpt(j,i) = a((i-1)*nzz+j,16)-tt(1,1);
                if rpt(j,i)<0
                    rpt(j,i) = 1e4;
                end
                if (abs(xx(j,i))<= L+2*ht && abs(zz(j,i)-zcenter)<= H +2*ht +hs )
                    ntag = ntag + 1;
                    rptout(ntag,1) = xx(j,i)*1e3;
                    rptout(ntag,2) = -(zz(j,i))*1e3;
                    rptout(ntag,3) = rpt(j,i);
                end
            else 
                rpt(j,i) = a((i-1)*nzz+j,col)/1e6;
            end
        end
    end


    figure(1)
    if col == 16
        contour(xx,zz,rpt,ll, 'ShowText' ,'on');colorbar;axis equal;hold on;
        delete(out);
        %%output the combined file
        fileID = fopen('a','w');
        fprintf(fileID,'%22.14e %22.14e %22.14e \n',rptout(:,:)');
        fclose(fileID);
        system(strcat('type'," ",'header\bp', num2str(bp),'\common_header.txt'," ",'>>'," ",out));
        system(strcat('type'," ",'header\bp', num2str(bp),'\rptheader.txt'," ",'>>'," ",out));
        system(strcat('type'," ",'a'," ",'>>'," ",out));
    else
        contourf(xx,zz,rpt);colorbar;axis equal;hold on;
    end
    delete tmp a;

end