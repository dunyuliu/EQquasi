function plotSlipStressProfileBP7(bp, mod, icend, mode)

    % The script creates required slip/stress_2/3_strike/depth.dat files from
    % EQquasi output required by the benchmark description of SEAS BP7.
    % 20221019.

    % mod             = 1; % model number
    % icend           = 3;
    % mode            = 1; % 1 strike; 2 depth profile
    % bp              = 7; % benchmark problem id.

    [path0,dx,nzz,l] = model_info_bp7(mod);

    % The outputs from EQquasi, p1/2output.txt contains 101 grids on each profile.
    % Write out the results from -400m to 400m, i.e. grid index 11 to 91. 
    if mode == 1
        x    = -500/1000:dx:500/1000; x=[0,0,x]; x = x*1000;
        o1   = 'slip_2_strike.dat';   o3 = 'slip_3_strike.dat';
        o2   = 'stress_2_strike.dat'; o4 = 'stress_3_strike.dat';
        file = 'p1output.txt';
        n    = floor(2*500/1000/dx) + 1; 
    elseif mode == 2
        x    = -500/1000:dx:500/1000; x=[0,0,x]; x = x*1000;
        o1   = 'slip_2_depth.dat';    o3 = 'slip_3_depth.dat';
        o2   = 'stress_2_depth.dat';  o4 = 'stress_3_depth.dat';
        file = 'p2output.txt';
        n    = floor(2*500/1000/dx) + 1;
    end

    o1 = strcat(path0, o1);
    o2 = strcat(path0, o2);
    o3 = strcat(path0, o3);
    o4 = strcat(path0, o4);

    tag       = 0; 
    tlast     = 0; 
    slip2last = zeros(1,n); 
    slip3last = zeros(1,n);
    for i = 1: icend
        path  = strcat(path0,'Q',num2str(i-1),'/',file);
        path2 = strcat(path0,'Q',num2str(i-1),'/global.dat');
        a     = load(path);
        b     = load(path2); nb=size(b,1);
        na    = size(a,1); 
        nstep = na/n;

        for k = 1:nstep
            t(k)           = a((k-1)*n+1,1);
            slip2(k,1:n)   = a((k-1)*n+1:k*n,2);
            slip3(k,1:n)   = a((k-1)*n+1:k*n,3);
            stress2(k,1:n) = a((k-1)*n+1:k*n,4)/1e6;
            stress3(k,1:n) = a((k-1)*n+1:k*n,5)/1e6;

            for m = 1 : nb
                if abs((t(k)-b(m,1)))<0.01
                    maxSlipRate(k) =log10(b(m,2));
                    break;
                end
            end
        end

        if i==1
            resslip2(:,1)       = t'; 
            resslip2(:,2)       = maxSlipRate';
            resslip2(:,3:3+n-1) = slip2(:,1:n);

            resslip3(:,1)       = t'; 
            resslip3(:,2)       = maxSlipRate';
            resslip3(:,3:3+n-1) = slip3(:,1:n);

            resstress2(:,1)     = t'; 
            resstress2(:,2)     = maxSlipRate';
            resstress2(:,3:3+n-1) = stress2(:,1:n);

            resstress3(:,1)     = t'; 
            resstress3(:,2)     = maxSlipRate';
            resstress3(:,3:3+n-1) = stress3(:,1:n);
        else
            resslip2tmp(:,1)    = t'+tlast; 
            resslip2tmp(:,2)    = maxSlipRate';
            for m = 1:nstep
                resslip2tmp(m,3:3+n-1)=slip2(m,1:n) + slip2last(1,1:n);
            end
            resslip2 = [resslip2;resslip2tmp];

            resslip3tmp(:,1)    = t'+tlast; 
            resslip3tmp(:,2)    = maxSlipRate';
            for m = 1:nstep
                resslip3tmp(m,3:3+n-1)=slip3(m,1:n) + slip3last(1,1:n);
            end
            resslip3 = [resslip3;resslip3tmp];

            resstress2tmp(:,1)  = t'+tlast; 
            resstress2tmp(:,2)  = maxSlipRate';
            resstress2tmp(:,3:3+n-1) = stress2(:,1:n); 
            resstress2          = [resstress2;resstress2tmp];      

            resstress3tmp(:,1)  = t'+tlast; 
            resstress3tmp(:,2)  = maxSlipRate';
            resstress3tmp(:,3:3+n-1) = stress3(:,1:n);
            resstress3          = [resstress3;resstress3tmp];    
        end
        tlast            = tlast + t(nstep);
        slip2last(1,1:n) = slip2last(1,1:n) + slip2(nstep,1:n);
        slip3last(1,1:n) = slip3last(1,1:n) + slip3(nstep,1:n);
        clear resslip2tmp a b slip2 slip3 stress2 stress3 maxSlipRate t;
        clear resslip3tmp resstress2tmp resstress3tmp;
    end

    ns2    = size(resslip2,1);
    ns3    = size(resslip3,1);
    nstr2  = size(resstress2,1);
    nstr3  = size(resstress3,1);

    figure(1)
    subplot(2,2,1)
    for i = 1 : ns2
        plot(resslip2(i,3:2+n));hold on;
    end

    subplot(2,2,2)
    for i = 1 : ns3
        plot(resslip3(i,3:2+n));hold on;
    end

    subplot(2,2,3)
    for i = 1 : nstr2
        plot(resstress2(i,3:2+n));hold on;
    end

    subplot(2,2,4)
    for i = 1 : nstr3
        plot(resstress3(i,3:2+n));hold on;
    end

    x1              = -400/1000:dx:400/1000; x1=[0,0,x1]; x1 = x1*1000;
    slip2(:,1:2)    = resslip2(:,1:2);
    slip2(:,3:83)   = resslip2(:,3+10:83+10);
    slip3(:,1:2)    = resslip3(:,1:2);
    slip3(:,3:83)   = resslip3(:,3+10:83+10);
    stress2(:,1:2)  = resstress2(:,1:2);
    stress2(:,3:83) = resstress2(:,3+10:83+10);
    stress3(:,1:2)  = resstress3(:,1:2);
    stress3(:,3:83) = resstress3(:,3+10:83+10);
    if mode == 1
        %%output the combined file
        delete(o1,o2,o3,o4);

        fileID = fopen('tmp','w');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',x1');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',slip2');
        fclose(fileID);
        system(strcat('type'," ",'header\bp', num2str(bp),'\common_header.txt'," ",'>>'," ",o1));
        system(strcat('type'," ",'header\bp', num2str(bp),'\slip2strikeheader.txt'," ",'>>'," ",o1));
        system(strcat('type'," ",'tmp'," ",'>>'," ",o1));

        fileID = fopen('tmp','w');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',x1');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',slip3');
        fclose(fileID);
        system(strcat('type'," ",'header\bp', num2str(bp),'\common_header.txt'," ",'>>'," ",o3));
        system(strcat('type'," ",'header\bp', num2str(bp),'\slip3strikeheader.txt'," ",'>>'," ",o3));
        system(strcat('type'," ",'tmp'," ",'>>'," ",o3));

        fileID = fopen('tmp','w');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',x1');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',stress2');
        fclose(fileID);
        system(strcat('type'," ",'header\bp', num2str(bp),'\common_header.txt'," ",'>>'," ",o2));
        system(strcat('type'," ",'header\bp', num2str(bp),'\stress2strikeheader.txt'," ",'>>'," ",o2));
        system(strcat('type'," ",'tmp'," ",'>>'," ",o2));

        fileID = fopen('tmp','w');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',x1');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',stress3');
        fclose(fileID);
        system(strcat('type'," ",'header\bp', num2str(bp),'\common_header.txt'," ",'>>'," ",o4));
        system(strcat('type'," ",'header\bp', num2str(bp),'\stress3strikeheader.txt'," ",'>>'," ",o4));
        system(strcat('type'," ",'tmp'," ",'>>'," ",o4));

    elseif mode == 2
        %%output the combined file
        delete(o1,o2,o3,o4);

        fileID = fopen('tmp','w');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',x1');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',slip2');
        fclose(fileID);
        system(strcat('type'," ",'header\bp', num2str(bp),'\common_header.txt'," ",'>>'," ",o1));
        system(strcat('type'," ",'header\bp', num2str(bp),'\slip2depthheader.txt'," ",'>>'," ",o1));
        system(strcat('type'," ",'tmp'," ",'>>'," ",o1));

        fileID = fopen('tmp','w');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',x1');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',slip3');
        fclose(fileID);
        system(strcat('type'," ",'header\bp', num2str(bp),'\common_header.txt'," ",'>>'," ",o3));
        system(strcat('type'," ",'header\bp', num2str(bp),'\slip3depthheader.txt'," ",'>>'," ",o3));
        system(strcat('type'," ",'tmp'," ",'>>'," ",o3));

        fileID = fopen('tmp','w');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',x1');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',stress2');
        fclose(fileID);
        system(strcat('type'," ",'header\bp', num2str(bp),'\common_header.txt'," ",'>>'," ",o2));
        system(strcat('type'," ",'header\bp', num2str(bp),'\stress2depthheader.txt'," ",'>>'," ",o2));
        system(strcat('type'," ",'tmp'," ",'>>'," ",o2));

        fileID = fopen('tmp','w');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',x1');
        fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',stress3');
        fclose(fileID);
        system(strcat('type'," ",'header\bp', num2str(bp),'\common_header.txt'," ",'>>'," ",o4));
        system(strcat('type'," ",'header\bp', num2str(bp),'\stress3depthheader.txt'," ",'>>'," ",o4));
        system(strcat('type'," ",'tmp'," ",'>>'," ",o4));    
    end
    delete tmp;
end