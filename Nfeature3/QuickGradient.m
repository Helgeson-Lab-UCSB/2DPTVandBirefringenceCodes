polyfitdeg=20;%degree of polynomial used for derivative fitting
gradbinsize=4;%gradient bin size (n x n square)
yres=1;
xres=yres;

%reinitialize gradient variables
dpux=zeros(imagexres*imageyres/(yres*yres),3);
dpvx=zeros(imagexres*imageyres/(yres*yres),3);
dpuy=zeros(imagexres*imageyres/(xres*xres),3);
dpvy=zeros(imagexres*imageyres/(xres*xres),3);
bindpux=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
bindpvx=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
bindpuy=zeros(imagexres*imageyres/(xres*xres*gradbinsize*gradbinsize),3);
bindpvy=zeros(imagexres*imageyres/(xres*xres*gradbinsize*gradbinsize),3);
FTP=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
MVGT=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
continuity=zeros(imagexres/(xres*gradbinsize),imageyres/(xres*gradbinsize));

tic
%begin gradient calculation 
disp(strcat('Calculating Gradients in x Direction...')); 

%derivatives of velocity (u and v) in x direction (fixed y's)
a=size(velocity);
o=1;
p=1;
for k=ceil(roiymin/yres):1:ceil(roiymax/yres)
    tempdx=[];
    for j=1:a(1)
        if  and(and(velocity(j,2)>=k*yres,velocity(j,2)<=(k+1)*yres),and(velocity(j,1)>=roixmin,velocity(j,1)<=roixmax))
            tempdx=[tempdx;velocity(j,:)];            
        end
    end
    tempdx=unique(tempdx,'rows'); %delete duplicate rows
    y=size(tempdx);
    if y(1)>polyfitdeg        
        ws = warning('off','all');  % Turn off warning
        [pux,Sux]=polyfit(tempdx(:,1),tempdx(:,3)*framerate,polyfitdeg); %fit polynomial to u
        [pvx,Svx]=polyfit(tempdx(:,1),tempdx(:,4)*framerate,polyfitdeg); %fit polynomial to v
        warning(ws);%turn warnings back on
        pux=polyder(pux);
        pvx=polyder(pvx);
        xprev=0;
        avgy=mean(tempdx(:,2));
        for x=transpose(tempdx(:,1))
            dpvx(p,1)=x;
            dpvx(p,2)=avgy;
            dpvx(p,3)=polyval(pvx,x);
            dpux(p,1)=x;
            dpux(p,2)=avgy;
            dpux(p,3)=polyval(pux,x);
            p=p+1;
        end
    end
end
%remove zeros from matrix
dpux(any(dpux==0,2),:)=[];
dpvx(any(dpvx==0,2),:)=[];
toc

disp(strcat('Calculating Gradients in y Direction...'));

%derivatives of velocity (u and v) in y direction (fixed x's)
a=size(velocity);
o=1;
p=1;
for k=ceil(roixmin/xres):1:ceil(roixmax/xres)
    tempdy=[];
    for j=1:a(1)
        if  and(and(velocity(j,1)>=k*xres,velocity(j,1)<=(k+1)*xres),and(velocity(j,2)>=roiymin,velocity(j,2)<=roiymax))
            tempdy=[tempdy;velocity(j,:)];            
        end
    end
    tempdy=unique(tempdy,'rows'); %delete duplicate rows
    y=size(tempdy);
    if y(1)>polyfitdeg
        ws = warning('off','all');  % Turn off warning
        puy=polyfit(tempdy(:,2),tempdy(:,3)*framerate,polyfitdeg); %fit polynomial to u
        pvy=polyfit(tempdy(:,2),tempdy(:,4)*framerate,polyfitdeg); %fit polynomial to v
        warning(ws);%turn warnings back on
        puy=polyder(puy);%take derivative of polynomials
        pvy=polyder(pvy);
        avgx=mean(tempdy(:,1));
        for y=transpose(tempdy(:,2))%loop over y's that contributed to fit
            dpvy(o,1)=avgx;
            dpvy(o,2)=y;
            dpvy(o,3)=polyval(pvy,y);
            o=o+1;
            dpuy(p,1)=avgx;
            dpuy(p,2)=y;
            dpuy(p,3)=polyval(puy,y);
            p=p+1;
        end
    end
end
toc

%remove zeros from matrix
dpuy(any(dpuy==0,2),:)=[];
dpvy(any(dpvy==0,2),:)=[];

disp(strcat('Binning du/dx...'));

%bin du/dx
a=size(dpux);
tic
matbindpux=zeros(imagexres/(yres*gradbinsize),imageyres/(yres*gradbinsize));
counter=zeros(imagexres/(yres*gradbinsize),imageyres/(yres*gradbinsize));
t=1;
for j=1:a(1)%pivot data to (x/binsize+1,y/binsize+1) matrix rounding to nearest integer and adding 1
    xclosest=round(dpux(j,1)/gradbinsize)+1;
    yclosest=round(dpux(j,2)/gradbinsize)+1;
    matbindpux(xclosest,yclosest)=matbindpux(xclosest,yclosest)+dpux(j,3);
    counter(xclosest,yclosest)=counter(xclosest,yclosest)+1;
end
for i=edgecutoff:ceil(imagexres/gradbinsize)-edgecutoff %average values in matrix elements
    for k=edgecutoff:ceil(imageyres/gradbinsize)-edgecutoff
        if counter(i+1,k+1)>0
            bindpux(t,1)=i*gradbinsize;
            bindpux(t,2)=k*gradbinsize;
            bindpux(t,3)=matbindpux(i+1,k+1)/counter(i+1,k+1);
            matbindpux(i+1,k+1)=matbindpux(i+1,k+1)/counter(i+1,k+1);
            t=t+1;
        end
    end
end
clear counter;
bindpux(any(bindpux==0,2),:)=[];%remove zeros from matrix
toc

disp(strcat('Binning dv/dx...'));

%bin dv/dx
a=size(dpvx);
tic
matbindpvx=zeros(imagexres/(yres*gradbinsize),imageyres/(yres*gradbinsize));
counter=zeros(imagexres/(yres*gradbinsize),imageyres/(yres*gradbinsize));
t=1;
for j=1:a(1)%pivot data to (x/binsize+1,y/binsize+1) matrix rounding to nearest integer and adding 1
    xclosest=round(dpvx(j,1)/gradbinsize)+1;
    yclosest=round(dpvx(j,2)/gradbinsize)+1;
    matbindpvx(xclosest,yclosest)=matbindpvx(xclosest,yclosest)+dpvx(j,3);
    counter(xclosest,yclosest)=counter(xclosest,yclosest)+1;
end
for i=edgecutoff:ceil(imagexres/gradbinsize)-edgecutoff %average values in matrix elements
    for k=edgecutoff:ceil(imageyres/gradbinsize)-edgecutoff
        if counter(i+1,k+1)>0
            bindpvx(t,1)=i*gradbinsize;
            bindpvx(t,2)=k*gradbinsize;
            bindpvx(t,3)=matbindpvx(i+1,k+1)/counter(i+1,k+1);
            matbindpvx(i+1,k+1)=matbindpvx(i+1,k+1)/counter(i+1,k+1);
            t=t+1;
        end
    end
end
clear counter;
bindpvx(any(bindpvx==0,2),:)=[];%remove zeros from matrix
toc

disp(strcat('Binning du/dy...'));

%bin du/dy
a=size(dpuy);
tic
matbindpuy=zeros(imagexres/(yres*gradbinsize),imageyres/(yres*gradbinsize));
counter=zeros(imagexres/(yres*gradbinsize),imageyres/(yres*gradbinsize));
t=1;
for j=1:a(1)%pivot data to (x/binsize+1,y/binsize+1) matrix rounding to nearest integer and adding 1
    xclosest=round(dpuy(j,1)/gradbinsize)+1;
    yclosest=round(dpuy(j,2)/gradbinsize)+1;
    matbindpuy(xclosest,yclosest)=matbindpuy(xclosest,yclosest)+dpuy(j,3);
    counter(xclosest,yclosest)=counter(xclosest,yclosest)+1;
end
for i=edgecutoff:ceil(imagexres/gradbinsize)-edgecutoff %average values in matrix elements
    for k=edgecutoff:ceil(imageyres/gradbinsize)-edgecutoff
        if counter(i+1,k+1)>0
            bindpuy(t,1)=i*gradbinsize;
            bindpuy(t,2)=k*gradbinsize;
            bindpuy(t,3)=matbindpuy(i+1,k+1)/counter(i+1,k+1);
            matbindpuy(i+1,k+1)=matbindpuy(i+1,k+1)/counter(i+1,k+1);
            t=t+1;
        end
    end
end
clear counter;
bindpuy(any(bindpuy==0,2),:)=[];%remove zeros from matrix
toc

disp(strcat('Binning dv/dy...'));

%bin dv/dy
a=size(dpvy);
tic
matbindpvy=zeros(imagexres/(yres*gradbinsize),imageyres/(yres*gradbinsize));
counter=zeros(imagexres/(yres*gradbinsize),imageyres/(yres*gradbinsize));
t=1;
for j=1:a(1)%pivot data to (x/binsize+1,y/binsize+1) matrix rounding to nearest integer and adding 1
    xclosest=round(dpvy(j,1)/gradbinsize)+1;
    yclosest=round(dpvy(j,2)/gradbinsize)+1;
    matbindpvy(xclosest,yclosest)=matbindpvy(xclosest,yclosest)+dpvy(j,3);
    counter(xclosest,yclosest)=counter(xclosest,yclosest)+1;
end
for i=edgecutoff:ceil(imagexres/gradbinsize)-edgecutoff %average values in matrix elements
    for k=edgecutoff:ceil(imageyres/gradbinsize)-edgecutoff
        if counter(i+1,k+1)>0
            bindpvy(t,1)=i*gradbinsize;
            bindpvy(t,2)=k*gradbinsize;
            bindpvy(t,3)=matbindpvy(i+1,k+1)/counter(i+1,k+1);
            matbindpvy(i+1,k+1)=matbindpvy(i+1,k+1)/counter(i+1,k+1);
            t=t+1;
        end
    end
end
clear counter;
bindpvy(any(bindpvy==0,2),:)=[];%remove zeros from matrix
toc

disp(strcat('Calculating FTP and VGTM...'));

%Flow type parameter (FTP) and magnitude of velocity gradient tensor (MVGT) calculation
MagE=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
MagO=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
FTP=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
MVGT=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
continuity=zeros(imagexres/(xres*gradbinsize),imageyres/(xres*gradbinsize));
v=1;
totalcont=0;
counter=0;
for i=edgecutoff:ceil(imagexres/gradbinsize)-edgecutoff
    for k=edgecutoff:ceil(imageyres/gradbinsize)-edgecutoff
        if matbindpux(i,k)~=0 && matbindpvx(i,k)~=0 && matbindpuy(i,k)~=0 && matbindpvy(i,k)~=0 %only evaluate when each tensor has a value at the point
            MagE(v,1)=i;
            MagE(v,2)=k;
            MagE(v,3)=sqrt(matbindpux(i,k).^2+matbindpvy(i,k).^2+0.5*(matbindpuy(i,k)+matbindpvx(i,k)).^2);
            MagO(v,1)=i;
            MagO(v,2)=k;
            MagO(v,3)=sqrt(0.25*((matbindpuy(i,k)-matbindpvx(i,k)).^2+(matbindpvx(i,k)-matbindpuy(i,k)).^2));
            FTP(v,1)=i;
            FTP(v,2)=k;
            FTP(v,3)=(MagE(v,3)-MagO(v,3))/(MagE(v,3)+MagO(v,3));
            MVGT(v,1)=i;
            MVGT(v,2)=k;
            MVGT(v,3)=sqrt(matbindpux(i,k).^2+matbindpuy(i,k).^2+matbindpvx(i,k).^2+matbindpvy(i,k).^2);
            continuity(v,1)=i;
            continuity(v,2)=k;
            continuity(v,3)=-(matbindpux(i,k)+matbindpvy(i,k));
            totalcont=totalcont+abs(continuity(v,3));
            counter=counter+1;
            v=v+1;
        end
    end
end
%remove zeros from matrix
FTP(any(FTP==0,2),:)=[];
MVGT(any(MVGT==0,2),:)=[];
toc

%graph FTP, MVGT, and extensional rate
figure
colormap jet;
scatter((FTP(:,1)*gradbinsize-xcenter)*lengthperpix,(FTP(:,2)*gradbinsize-ycenter)*lengthperpix,gradbinsize*4,FTP(:,3),'filled','s')
caxis([-1,1])
colorbar;axis equal;axis([(roixmin-xcenter)*lengthperpix (roixmax-xcenter)*lengthperpix (roiymin-ycenter)*lengthperpix (roiymax-ycenter)*lengthperpix]);
xlabel(lengthunit)
ylabel(lengthunit)
title(colorbar,strcat('FTP'))

figure
colormap jet;
scatter((MVGT(:,1)*gradbinsize-xcenter)*lengthperpix,(MVGT(:,2)*gradbinsize-ycenter)*lengthperpix,gradbinsize*4,MVGT(:,3),'filled','s')
caxis([0,maxMVGT])
colorbar;axis equal;axis([(roixmin-xcenter)*lengthperpix (roixmax-xcenter)*lengthperpix (roiymin-ycenter)*lengthperpix (roiymax-ycenter)*lengthperpix]);
xlabel(lengthunit)
ylabel(lengthunit)
title(colorbar,strcat('MVGT (1/',timeunit,')'))

figure
colormap jet;
scatter((MagE(:,1)*gradbinsize-xcenter)*lengthperpix,(MagE(:,2)*gradbinsize-ycenter)*lengthperpix,gradbinsize*4,MagE(:,3),'filled','s')
caxis([-1,1])
colorbar;axis equal;axis([(roixmin-xcenter)*lengthperpix (roixmax-xcenter)*lengthperpix (roiymin-ycenter)*lengthperpix (roiymax-ycenter)*lengthperpix]);
xlabel(lengthunit)
ylabel(lengthunit)
title(colorbar,strcat('||E|| (1/',timeunit,')'))

% figure
% colormap jet;
% scatter((bindpvy(:,1)-xcenter)*lengthperpix,(bindpvy(:,2)-ycenter)*lengthperpix,gradbinsize*2,bindpvy(:,3),'filled','s')
% caxis([0,maxMVGT])
% colorbar;axis equal;axis([(roixmin-xcenter)*lengthperpix (roixmax-xcenter)*lengthperpix (roiymin-ycenter)*lengthperpix (roiymax-ycenter)*lengthperpix]);
% xlabel(lengthunit)
% ylabel(lengthunit)
% title(colorbar,strcat('Extensional Rate (1/',timeunit,')'))
% 
% figure
% colormap jet;
% scatter((bindpux(:,1)-xcenter)*lengthperpix,(bindpux(:,2)-ycenter)*lengthperpix,gradbinsize*4,bindpux(:,3),'filled','s')
% caxis([0,maxMVGT])
% colorbar;axis equal;axis([(roixmin-xcenter)*lengthperpix (roixmax-xcenter)*lengthperpix (roiymin-ycenter)*lengthperpix (roiymax-ycenter)*lengthperpix]);
% xlabel(lengthunit)
% ylabel(lengthunit)
% title(colorbar,strcat('Extensional Rate (1/',timeunit,')'))

% figure
% colormap jet;
% scatter((continuity(:,1)*gradbinsize-xcenter)*lengthperpix,(continuity(:,2)*gradbinsize-ycenter)*lengthperpix,gradbinsize*4,continuity(:,3),'filled','s')
% caxis([-3,3])
% colorbar;axis equal;axis([(roixmin-xcenter)*lengthperpix (roixmax-xcenter)*lengthperpix (roiymin-ycenter)*lengthperpix (roiymax-ycenter)*lengthperpix]);
% xlabel(lengthunit)
% ylabel(lengthunit)
% title(colorbar,strcat('dw/dz (1/',timeunit,')'))


%save(strcat(imageloc,'\',imagename,'matlabdata.mat'))%save generated matlab data