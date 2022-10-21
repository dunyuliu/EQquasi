function plotGlobaldat(bp, mod, icend)
% bp     = 7;
% mod    = 1;
% icend  = 3;
    txt = 'Generating global.dat ... ...'
    if     bp == 5 
        [path,dx,nzz,l] = model_info_bp5(mod);
    elseif bp == 7
        [path,dx,nzz,l] = model_info_bp7(mod);
    end

    out      = strcat(path,'global.dat');

    finalt    = 0; 
    finald    = 0;
    for ic = 1:icend
        path0 = strcat(path,'Q',num2str(ic-1),'/');
        glo   = load(strcat(path0,'global.dat'));
        % 1-time 2-maxsliprate 3-momrate 4-tao*A 5-slip*A 6-A
        glo(:,1) = glo(:,1) + finalt;
        n     = size(glo,1);

        if ic ==1
            totglo = glo;
            finalt = glo(n,1);
        elseif ic >1
            finalt = glo(n,1);
            totglo =[totglo;glo;];
        end

    end
    totglo(:,2) = log10(totglo(:,2));

    h1 = figure(1);
    set(h1,'position',[100 100 700 500]);
    subplot(4,1,1)
    plot(totglo(:,1),totglo(:,2)); title('Max Sliprate (m/s)');
    subplot(4,1,2)
    plot(totglo(:,1),totglo(:,3)); title('Moment Rate (N-m/s)');
    subplot(4,1,3)
    plot(totglo(:,1),totglo(:,4)); title('Ave Shear (MPa)');
    subplot(4,1,4)
    plot(totglo(:,1),totglo(:,5)); title('Ave Slip (m)');

    %%output the combined file
    delete(out)
    fileID = fopen('tmp','w');
    if bp == 5
        fprintf(fileID,'%22.14e %15.7e %15.7e\n',totglo(:,1:3)');
    elseif bp == 7
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e\n',totglo(:,1:4)');
    end

    fclose(fileID);
    system(strcat('type'," ",'header\bp', num2str(bp),'\common_header.txt'," ",'>>'," ",out));
    system(strcat('type'," ",'header\bp', num2str(bp),'\sourceheader.txt'," ",'>>'," ",out));
    system(strcat('type'," ",'tmp'," ",'>>'," ",out));
    delete tmp;
end