clear all;close all;


dx = 1;
%M = dlmread('slip.txt',' ',[(t1-1)*n 0 (t2-1)*n-1 3]);
%load('MyColormaps.mat','rgb');
ModelName=strcat('SEAS-BP4-M1');
for ic = 1:1
    FolderName=strcat('./C',num2str(ic-1),'/');    
    M=load(strcat(FolderName,'slip.txt'));
    dataq=load(strcat(FolderName,'global.dat'));
    tt=load(strcat(FolderName,'tdyna.txt'));
    tmp = dataq(:,1); 
    ntmp = size(dataq,1);
    tag = 0;
    for i = 1:ntmp
        if (tmp(i,1)<=tt(1,1))
            if mod(i,20) == 1
                tag = tag + 1;
                T(tag,1) = tmp(i,1);
            end
        else
            if mod(i,10) == 1
                tag = tag + 1;
                T(tag,1) = tmp(i,1);
            end
        end
    end
    T=T/60/60/24/365;
    position=[10 10 1200 600];
    set(0,'DefaultFigurePosition', position);
    fig1=figure(1);
    winsize= get(fig1,'Position');
    winsize(1:2)=[0 0];
    %A=moviein(nstep,fig1,winsize);
    set(fig1,'NextPlot','replacechildren');
    %fid = fopen('slip.txt');
    x=-60:1:60;z=-120:1:0;
    nx=size(x',1);nz=size(z',1);
    [xx,zz]=meshgrid(x,z);
    n=nx*nz;
    [totstep,tmp1]=size(T);

    a0=T(:,1);st=zeros(totstep,3);tn=zeros(totstep,3);ts=zeros(totstep,3);
    st1=[0, -60];st2=[-20, -70];st3=[20,-70];
    label=0;
    for nt=1:totstep %+floor(3/4*totstep):totstep
    %         if label==1
    %             break
    %         end
        nt
        t1=nt;
        t2=t1+1;
        for i=1:nx
            for j=1:nz
                sl(j,i)=M((t1-1)*n+(i-1)*nz+j,1);
                sr(j,i)=M((t1-1)*n+(i-1)*nz+j,2);
                if sr(j,i)<0
                    sr(j,i)=1000;
                    label=1;
                    xcoor=(i-1)*dx-60
                    zcoor=(j-1)*dx-120
                end
                sr(j,i)=log10(sr(j,i));
                psi(j,i)=log10(M((t1-1)*n+(i-1)*nz+j,4));
                tao(j,i)=M((t1-1)*n+(i-1)*nz+j,5);
                tnrm(j,i)=M((t1-1)*n+(i-1)*nz+j,6);
    %                 if (tao(j,i)<10e6)
    %                     tao(j,i)=30e6;
    %                 end
                    if (abs((i-1)*dx-st1(1)-60)<0.001 &&abs((j-1)*dx-120-st1(2))<0.001 )
                    st(nt,1)=sr(j,i);
                    ts(nt,1)=tao(j,i)/1d6;
                    tn(nt,1)=tnrm(j,i)/1d6;   
                    elseif (abs((i-1)*dx-st2(1)-60)<0.001 &&abs((j-1)*dx-120-st2(2))<0.001 )
                    st(nt,2)=sr(j,i);
                    ts(nt,2)=tao(j,i)/1d6;
                    tn(nt,2)=tnrm(j,i)/1d6;   
                    elseif (abs((i-1)*dx-st3(1)-60)<0.001 &&abs((j-1)*dx-120-st3(2))<0.001 )
                    st(nt,3)=sr(j,i); 
                    ts(nt,3)=tao(j,i)/1d6;
                    tn(nt,3)=tnrm(j,i)/1d6;   
                end
            end
        end
        [row, col] = find(isnan(sr))
        figure (fig1)

        subplot('position',[0.05 0.5 0.28 0.39]);
        contourf(xx,zz,tnrm/1e6,'linestyle','none');colorbar;colormap(jet);caxis([-80,-40]);
        axis equal;
        title(strcat(ModelName,' - tnrm (MPa) -',num2str(T(nt,1)),' yr'));

        subplot('position',[0.05+0.3 0.5 0.28 0.39]);
        contourf(xx,zz,sr,'linestyle','none');title('sliprate (m/s)');
        colorbar;colormap(jet);caxis([-14,2]);
        hold on; plot(st1(1),st1(2),'b*');
        hold on; plot(st2(1),st2(2),'k*');
        hold on; plot(st3(1),st3(2),'m*');
        axis equal;

        subplot('position',[0.05 0.05 0.28 0.39]);
        contourf(xx,zz,tao/1e6,'linestyle','none');title('tstk (MPa)');colorbar;colormap(jet);caxis([20,60]);%colorbar;xlim([-15,15]);ylim([-18,0]);
        axis equal;

        subplot('position',[0.05+0.3 0.05 0.28 0.39]);
        contourf(xx,zz,psi,'linestyle','none');title('psi');colorbar;colormap(jet);caxis([1,13]);
        axis equal;

       subplot('position',[0.05+0.63 0.1 0.28 0.22]);
        plot(a0(1:nt),st(1:nt,1),'b');hold on;
        plot(a0(1:nt),st(1:nt,2),'k');hold on;
        plot(a0(1:nt),st(1:nt,3),'m');hold on;
        plot(a0(nt,1),st(nt,1),'b*');hold on;
        plot(a0(nt,1),st(nt,2),'k*');hold on;
        plot(a0(nt,1),st(nt,3),'m*');
        ylim([-16 2]);xlim([0 T(totstep,1)]);
        ylabel('Log(sliprate)(log m/s) at STs');xlabel('Years');

                subplot('position',[0.05+0.63 0.4 0.28 0.22]);
        plot(a0(1:nt),ts(1:nt,1),'b');hold on;
        plot(a0(1:nt),ts(1:nt,2),'k');hold on;
        plot(a0(1:nt),ts(1:nt,3),'m');hold on;
        plot(a0(nt,1),ts(nt,1),'b*');hold on;
        plot(a0(nt,1),ts(nt,2),'k*');hold on;
        plot(a0(nt,1),ts(nt,3),'m*');
        ylim([20 60]);xlim([0 T(totstep,1)]);
        ylabel('tstk (MPa) at STs');

         subplot('position',[0.05+0.63 0.7 0.28 0.22]);
        plot(a0(1:nt),tn(1:nt,1),'b');hold on;
        plot(a0(1:nt),tn(1:nt,2),'k');hold on;
        plot(a0(1:nt),tn(1:nt,3),'m');hold on;
        plot(a0(nt,1),tn(nt,1),'b*');hold on;
        plot(a0(nt,1),tn(nt,2),'k*');hold on;
        plot(a0(nt,1),tn(nt,3),'m*');
        ylim([-80 -40]);xlim([0 T(totstep,1)]);
        ylabel('tnrm (MPa) at STs'); legend('x=0 km','x=-20 km','x=20 km');      

        A(nt)=getframe(fig1,winsize); 

    end   
end
save(strcat('./0Movie/',ModelName,'.mat'),'A','-v7.3');

