clear all; close all;

out = 'global.txt';
finalt = 0;
for ic = 1:1
    path = strcat('../work/C',num2str(ic-1),'/')
    glo = load(strcat(path,'global.dat'));
    glo(:,1) = glo(:,1) + finalt;
    n = size(glo,1);

    if ic ==1
        totglo = glo;
        finalt = glo(n,1);
    elseif ic >1
        finalt = glo(n,1);
        totglo =[totglo;glo;];
    end
        
end
totglo(:,2) = log10(totglo(:,2));

np = size(totglo,1);
for i = 1:np
    if totglo(i,6) > 0 
        totglo(i,4) =  totglo(i,4)/ totglo(i,6)/1e6;
        totglo(i,5) =  totglo(i,5)/ totglo(i,6);
    else
        totglo(i,4) = 0;
        totglo(i,5) = 0;
    end
end

h1 = figure(1);
set(h1,'position',[100 100 700 500]);
subplot(4,1,1)
plot(totglo(:,1),totglo(:,2)); title('Max Sliprate s (m/s)');
subplot(4,1,2)
plot(totglo(:,1),totglo(:,3)); title('Moment rate (N-m/s)');
subplot(4,1,3)
plot(totglo(:,1),totglo(:,4)); title('Ave Shear MPa');
subplot(4,1,4)
plot(totglo(:,1),totglo(:,5)); title('Ave slip m');

%%output the combined file
system(strcat('rm'," ",out));
fileID = fopen('tmp','w');
fprintf(fileID,'%22.14e %15.7e %15.7e\n',totglo(:,1:3)');
fclose(fileID);
system(strcat('touch'," ",out));
system(strcat('cat'," ",'sourceheader.txt'," ",'>>'," ",out));
system(strcat('cat'," ",'tmp'," ",'>>'," ",out));