clear all; close all;
% TIME(1);SLIPS(2);SLIPD(3);SLIPRS(4),SLIPRD(5);TSTK(6);TDIP(7);STATE(8)
mod = 11;
icend = 4;
col = 'k';

[path,dx,nzz,l] = model_info(mod);
path

for s = 1:10
if s==1
    st1 = 'fltst_strk000dp000.txt'; outst1 = 'fltst_strk+00dp+00.txt';
elseif s == 2
    st1 = 'fltst_strk000dp010.txt'; outst1 = 'fltst_strk+00dp+10.txt';
elseif s == 3
    st1 = 'fltst_strk000dp022.txt'; outst1 = 'fltst_strk+00dp+22.txt';    
elseif s == 4
    st1 = 'fltst_strk016dp000.txt'; outst1 = 'fltst_strk+16dp+00.txt';
elseif s == 5
    st1 = 'fltst_strk016dp010.txt'; outst1 = 'fltst_strk+16dp+10.txt';
elseif s == 6
    st1 = 'fltst_strk036dp000.txt'; outst1 = 'fltst_strk+36dp+00.txt';    
elseif s == 7
    st1 = 'fltst_strk-016dp000.txt'; outst1 = 'fltst_strk-16dp+00.txt';    
elseif s == 8
    st1 = 'fltst_strk-016dp010.txt'; outst1 = 'fltst_strk-16dp+10.txt';  
elseif s == 9
    st1 = 'fltst_strk-024dp010.txt'; outst1 = 'fltst_strk-24dp+10.txt';    
elseif s == 10
    st1 = 'fltst_strk-036dp000.txt'; outst1 = 'fltst_strk-36dp+00.txt'; 
end
%st1 = 'fltst_strk--220dp070.txt'; outst1 = 'fltst_strk-225dp+750.txt';
out = strcat(path, 'res/', outst1);

finalt = 0;
slips = 0; slipd = 0;
for ic = 1 : icend
    path1 = strcat(path,'Q',num2str(ic-1),'/');
    glo = load(strcat(path1,'global.dat'));
    tt = load(strcat(path1,'tdyna.txt'));
    res = load(strcat(path1,st1));
    res(:,1) = res(:,1) + finalt;
    res(:,2) = res(:,2) + slips;
    res(:,3) = res(:,3) + slipd;
    n = size(res,1);
%     n = size(glo,1);
%     for i = 1: n
%         if glo(i,1)==tt(1,1)
%             ntagstart = i;
%         end
%         if glo(i,1) ==tt(1,2)
%             ntagend = i;
%         end
%     end
    if ic ==1
        totres = res;
        finalt = res(n,1);
        slips = res(n,2);
        slipd = res(n,3);
    elseif ic >1
        finalt = res(n,1);
        slips = res(n,2);
        slipd = res(n,3);        
        totres =[totres;res;];
    end
        
end
totres1(:,1:3) = totres(:,1:3);
totres1(:,4:5) = log10(totres(:,4:5));
totres1(:,6:7) = totres(:,6:7);
totres1(:,8) = totres(:,8);
h1 = figure(1);
set(h1,'position',[100 100 700 500]);
subplot(4,1,1)
plot(totres(:,1),log10(totres(:,4)),col); title('Sliprate s (m/s)');hold on;
subplot(4,1,2)
plot(totres(:,1),log10(abs(totres(:,5))),col); title('Sliprate d (m/s)');hold on;
subplot(4,1,3)
plot(totres(:,1),totres(:,6),col); title('Shear Stress (MPa)');hold on;
% subplot(4,1,4)
% plot(totres(:,1),totres(:,7),col); title('Normal Stress (MPa)');hold on;


%%output the combined file
delete(out);
fileID = fopen('tmp','w');
fprintf(fileID,'%22.14e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e\n',totres1(:,:)');
fclose(fileID);
system(strcat('type'," ",'.\header\stheader.txt'," ",'>'," ",out));
system(strcat('type'," ",'tmp'," ",'>>'," ",out));
% system(strcat('touch'," ",out));
% system(strcat('cat'," ",'stheader.txt'," ",'>>'," ",out));
% system(strcat('cat'," ",'tmp'," ",'>>'," ",out));

clear glo tt res totres;
end