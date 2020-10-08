clear all; close all;

for ic = 1:1
    path = strcat('../work/C',num2str(ic-1),'/')
    glo = load(strcat(path,'global.dat'));
    tt = load(strcat(path,'tdyna.txt'));
    n = size(glo,1);
    for i = 1: n
        if glo(i,1)==tt(1,1)
            ntagstart = i;
        end
        if glo(i,1) ==tt(1,2)
            ntagend = i;
        end
    end
    if ic ==1 
        res = glo;
    else
        res = [res;glo];
    end
end

h1 = figure(1);
set(h1,'position',[100 100 700 500]);
subplot(4,1,1)
plot(res(:,1),log10(res(:,2))); title('Maxsliprate (m/s)');
subplot(4,1,2)
plot(res(:,1),res(:,3)); title('Momment rate (N-m/s)');
subplot(4,1,3)
plot(res(ntagstart:ntagend,1)-tt(1,1),res(ntagstart:ntagend,4)/1e6); 
ylim([20 40]);
title('Average shear stress (MPa)');
subplot(4,1,4)
plot(res(ntagstart:ntagend,1)-tt(1,1),res(ntagstart:ntagend,5)); 
title('Average slip (m)');