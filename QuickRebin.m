gradbinsize=16;

disp(strcat('Binning du/dx...'));

%bin du/dx
a=size(dpux);
tic
matbindpux=zeros(imagexres/(gradbinsize)+1,imageyres/(gradbinsize)+1);
bindpux=zeros(imagexres*imageyres/(gradbinsize*gradbinsize),3);
counter=zeros(imagexres/(gradbinsize)+1,imageyres/(gradbinsize)+1);
t=1;
for j=1:a(1)%pivot data to (x/binsize+1,y/binsize+1) matrix rounding to nearest integer and adding 1
    xclosest=round(dpux(j,1)/gradbinsize)+1;
    yclosest=round(dpux(j,2)/gradbinsize)+1;
    matbindpux(xclosest,yclosest)=matbindpux(xclosest,yclosest)+dpux(j,3)*dpux(j,4);
    counter(xclosest,yclosest)=counter(xclosest,yclosest)+dpux(j,4);
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
matbindpvx=zeros(imagexres/(gradbinsize)+1,imageyres/(gradbinsize)+1);
bindpvx=zeros(imagexres*imageyres/(gradbinsize*gradbinsize),3);
counter=zeros(imagexres/(gradbinsize)+1,imageyres/(gradbinsize)+1);
t=1;
for j=1:a(1)%pivot data to (x/binsize+1,y/binsize+1) matrix rounding to nearest integer and adding 1
    xclosest=round(dpvx(j,1)/gradbinsize)+1;
    yclosest=round(dpvx(j,2)/gradbinsize)+1;
    matbindpvx(xclosest,yclosest)=matbindpvx(xclosest,yclosest)+dpvx(j,3)*dpvx(j,4);
    counter(xclosest,yclosest)=counter(xclosest,yclosest)+dpvx(j,4);
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
matbindpuy=zeros(imagexres/(gradbinsize)+1,imageyres/(gradbinsize)+1);
bindpuy=zeros(imagexres*imageyres/(gradbinsize*gradbinsize),3);
counter=zeros(imagexres/(gradbinsize)+1,imageyres/(gradbinsize)+1);
t=1;
for j=1:a(1)%pivot data to (x/binsize+1,y/binsize+1) matrix rounding to nearest integer and adding 1
    xclosest=round(dpuy(j,1)/gradbinsize)+1;
    yclosest=round(dpuy(j,2)/gradbinsize)+1;
    matbindpuy(xclosest,yclosest)=matbindpuy(xclosest,yclosest)+dpuy(j,3)*dpuy(j,4);
    counter(xclosest,yclosest)=counter(xclosest,yclosest)+dpuy(j,4);
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
matbindpvy=zeros(imagexres/(gradbinsize)+1,imageyres/(gradbinsize)+1);
bindpvy=zeros(imagexres*imageyres/(gradbinsize*gradbinsize),3);
counter=zeros(imagexres/(gradbinsize)+1,imageyres/(gradbinsize)+1);
t=1;
for j=1:a(1)%pivot data to (x/binsize+1,y/binsize+1) matrix rounding to nearest integer and adding 1
    xclosest=round(dpvy(j,1)/gradbinsize)+1;
    yclosest=round(dpvy(j,2)/gradbinsize)+1;
    matbindpvy(xclosest,yclosest)=matbindpvy(xclosest,yclosest)+dpvy(j,3)*dpvy(j,4);
    counter(xclosest,yclosest)=counter(xclosest,yclosest)+dpvy(j,4);
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
MagE=zeros(imagexres*imageyres/(gradbinsize*gradbinsize),3);
MagO=zeros(imagexres*imageyres/(gradbinsize*gradbinsize),3);
FTP=zeros(imagexres*imageyres/(gradbinsize*gradbinsize),3);
MVGT=zeros(imagexres*imageyres/(gradbinsize*gradbinsize),3);
continuity=zeros(imagexres/(gradbinsize)+1,imageyres/(gradbinsize)+1);
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
scatter((FTP(:,1)*gradbinsize-xcenter)*lengthperpix,(FTP(:,2)*gradbinsize-ycenter)*lengthperpix,gradbinsize,FTP(:,3),'filled','s')
caxis([-1,1])
colorbar;axis equal;axis([(0-xcenter)*lengthperpix (imagexres-xcenter)*lengthperpix (0-ycenter)*lengthperpix (imageyres-ycenter)*lengthperpix]);
xlabel(lengthunit)
ylabel(lengthunit)
title(colorbar,strcat('FTP'))

figure
colormap jet;
scatter((MVGT(:,1)*gradbinsize-xcenter)*lengthperpix,(MVGT(:,2)*gradbinsize-ycenter)*lengthperpix,gradbinsize,MVGT(:,3),'filled','s')
caxis([0,maxMVGT])
colorbar;axis equal;axis([(0-xcenter)*lengthperpix (imagexres-xcenter)*lengthperpix (0-ycenter)*lengthperpix (imageyres-ycenter)*lengthperpix]);
xlabel(lengthunit)
ylabel(lengthunit)
title(colorbar,strcat('MVGT (1/',timeunit,')'))

% figure
% colormap jet;
% scatter((bindpvy(:,1)-xcenter)*lengthperpix,(bindpvy(:,2)-ycenter)*lengthperpix,gradbinsize,bindpvy(:,3),'filled','s')
% caxis([0,maxMVGT])
% colorbar;axis equal;axis([(0-xcenter)*lengthperpix (imagexres-xcenter)*lengthperpix (0-ycenter)*lengthperpix (imageyres-ycenter)*lengthperpix]);
% xlabel(lengthunit)
% ylabel(lengthunit)
% title(colorbar,strcat('Extensional Rate (1/',timeunit,')'))
% 
% figure
% colormap jet;
% scatter((bindpux(:,1)-xcenter)*lengthperpix,(bindpux(:,2)-ycenter)*lengthperpix,gradbinsize*4,bindpux(:,3),'filled','s')
% caxis([0,maxMVGT])
% colorbar;axis equal;axis([(0-ycenter)*lengthperpix (imageyres-ycenter)*lengthperpix (0-ycenter)*lengthperpix (imageyres-ycenter)*lengthperpix]);
% xlabel(lengthunit)
% ylabel(lengthunit)
% title(colorbar,strcat('Extensional Rate (1/',timeunit,')'))
% 
% figure
% colormap jet;
% scatter((continuity(:,1)*gradbinsize-xcenter)*lengthperpix,(continuity(:,2)*gradbinsize-ycenter)*lengthperpix,gradbinsize*4,continuity(:,3),'filled','s')
% caxis([-3,3])
% colorbar;axis equal;axis([(0-ycenter)*lengthperpix (imageyres-ycenter)*lengthperpix (0-ycenter)*lengthperpix (imageyres-ycenter)*lengthperpix]);
% xlabel(lengthunit)
% ylabel(lengthunit)
% title(colorbar,strcat('dw/dz (1/',timeunit,')'))

save(strcat(imageloc,'\',imagename,'matlabdata.mat'))%save generated matlab data