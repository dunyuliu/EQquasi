function plotOffFaultSt(bp, mod, icend)
    txt = 'Generating OffFaultStation.txt ... ...'
    % TIME(1);SLIPS(2);SLIPD(3);SLIPRS(4),SLIPRD(5);TSTK(6);TDIP(7);STATE(8)
    % bp     = 7;
    % mod    = 1;
    % icend  = 3;

    col    = 'k';

    if     bp == 5 
        [path,dx,nzz,l] = model_info_bp5(mod);
        NumTotSt = 9;
    elseif bp == 7
        [path,dx,nzz,l] = model_info_bp7(mod);
        NumTotSt = 5;
    end

    for s = 1:NumTotSt
        if bp == 5
            if s==1
                st1 = 'srfst_strk000st008dp000.txt'; outst1 = 'blkst_strk+00fn+08dp+00.txt';
            elseif s == 2
                st1 = 'srfst_strk000st008dp010.txt'; outst1 = 'blkst_strk+00dp+08dp+10.txt';
            elseif s == 3
                st1 = 'srfst_strk000st016dp000.txt'; outst1 = 'blkst_strk+00fn+16dp+00.txt';    
            elseif s == 4
                st1 = 'srfst_strk000st016dp010.txt'; outst1 = 'blkst_strk+00fn+16dp+10.txt';
            elseif s == 5
                st1 = 'srfst_strk000st032dp000.txt'; outst1 = 'blkst_strk+00fn+32dp+00.txt';
            elseif s == 6
                st1 = 'srfst_strk000st032dp010.txt'; outst1 = 'blkst_strk+00fn+32dp+10.txt';    
            elseif s == 7
                st1 = 'srfst_strk000st048dp000.txt'; outst1 = 'blkst_strk00fn+48dp+00.txt';  
            elseif s == 8
                st1 = 'srfst_strk016st008dp000.txt'; outst1 = 'blkst_strk+16fn+08dp+00.txt';    
            elseif s == 9
                st1 = 'srfst_strk-016st008dp000.txt'; outst1 = 'blkst_strk-16fn+08dp+00.txt';        
            end
        elseif bp == 7
            if s==1
                st1 = 'srfst_strk000st200dp000.txt'; outst1 = 'blkst_strk+000fn+200dp+000.txt';
            elseif s == 2
                st1 = 'srfst_strk000st400dp000.txt'; outst1 = 'blkst_strk+000dp+400dp+000.txt';
            elseif s == 3
                st1 = 'srfst_strk000st400dp300.txt'; outst1 = 'blkst_strk+000fn+400dp+300.txt';    
            elseif s == 4
                st1 = 'srfst_strk300st400dp000.txt'; outst1 = 'blkst_strk+300fn+400dp+000.txt';
            elseif s == 5
                st1 = 'srfst_strk-300st400dp000.txt'; outst1 = 'blkst_strk-300fn+400dp+000.txt';
            end
        end    

    out = strcat(path, outst1);

    finalt = 0;
    ds = 0; dd = 0; dfn = 0;
    for ic = 1 : icend
        path1 = strcat(path,'Q',num2str(ic-1),'/');
        glo   = load(strcat(path1,'global.dat'));
        tt    = load(strcat(path1,'tdyna.txt'));
        res   = load(strcat(path1,st1));
        res(:,1) = res(:,1) + finalt;   %time
        res(:,2) = res(:,2) + ds;       %h-disp positive x+
        res(:,4) = res(:,4) + dd;       %v-disp positive z-
        res(:,6) = res(:,6) + dfn;      %fn-disp positive y+

        n = size(res,1);
        
        if ic ==1
            totres = res;
            finalt = res(n,1);
            ds = res(n,2);
            dd = res(n,4);
            dfn = res(n,6);
        elseif ic >1
            finalt = res(n,1);
            ds = res(n,2);
            dd = res(n,4);
            dfn = res(n,6);    
            totres =[totres;res;];
        end

    end
    % Quants in final output
    % 1-fnd 2-sd 3-vd(positive downward) 4=fnvel 5-svel 6-vvel
    totres1(:,1) = totres(:,1);
    totres1(:,2) = totres(:,6);
    totres1(:,3) = totres(:,2);
    totres1(:,4) = totres(:,4);

    % Calculate particle velocity from particle disps.
    nsample = size(totres1,1);
    totres1(1,5:7) = 0;
    for i = 2:nsample
        dt = totres1(i,1) - totres1(i-1,1);
        totres1(i,5) = (totres1(i,2) - totres1(i-1,2))/dt;
        totres1(i,6) = (totres1(i,3) - totres1(i-1,3))/dt;
        totres1(i,7) = (totres1(i,4) - totres1(i-1,4))/dt;
    end

    % totres1(:,5) = log10(abs(totres(:,7)+1e-30)); 
    % totres1(:,6) = log10(abs(totres(:,3)+1e-30));
    % totres1(:,7) = log10(abs(totres(:,5)+1e-30));

    totres1(1,2:7) = 1e-30;

    h1 = figure(1);
    set(h1,'position',[100 100 700 500]);
    subplot(6,1,1)
    plot(totres1(:,1),totres1(:,2),col); title('disp fn (m)');hold on;
    subplot(6,1,2)
    plot(totres1(:,1),totres1(:,3),col); title('disp str (m)');hold on;
    subplot(6,1,3)
    plot(totres1(:,1),totres1(:,4),col); title('disp dip (m)');hold on;

    subplot(6,1,4)
    plot(totres1(:,1),totres1(:,5),col); title('log10 vel fn (m/s)');hold on;
    subplot(6,1,5)
    plot(totres1(:,1),totres1(:,6),col); title('log10 vel str (m/s)');hold on;
    subplot(6,1,6)
    plot(totres1(:,1),totres1(:,7),col); title('log10 vel dip');hold on;

    % %%output the combined file
    delete(out);
    fileID = fopen('tmp','w');
    fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n',totres1(:,:)');
    fclose(fileID);
    system(strcat('type'," ",'header\bp',num2str(bp),'\common_header.txt'," ",'>>'," ",out));
    system(strcat('type'," ",'header\bp',num2str(bp),'\offstheader.txt'," ",'>>'," ",out));
    system(strcat('type'," ",'tmp'," ",'>>'," ",out));

    clear glo tt res totres;
    end
    delete tmp;
end