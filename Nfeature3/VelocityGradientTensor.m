%Analyze set of processed images for velocity, FTP, and MVGT fields
%functions required:
%Nfeature3.m (and functions that this one uses)
%trackmem.m (and functions that this one uses)

%inputs---------------------------------------------

%image parameters
lengthperpix=2/384;%distance per pixel
lengthunit='mm'; %unit of distance
framerate=15;%frames per time
timeunit='sec'; %unit of time
imagexres=2048;%pixel resolution of images
imageyres=1088;
numbofim=1000;%number of images
imageloc='C:\Users\Patrick\16-6-10\NIST WLM H=5mm W=2mm XSL\0.24mlmin15fps';%image files location
imagename='0.24mlmin(2)';%image name without '0000','0001',etc or '.tif'

%graphing parameters
xcenter=912;%pixel center of geometry
ycenter=568;
maxMVGT=0.4;%colorbar maximum for MVGT (1/time)

%velocity and gradient quality parameters
trajectorygap=2; %images between velocity calculation
yres=1;%y slice in derivative calculation (please leave as 1)
xres=yres; %x slice in derivative calculation (please leave as 1)
polyfitdeg=16;%degree of polynomial used for derivative fitting
roixmin=580;
roixmax=1250;
roiymin=250;
roiymax=950;

%binning parameters
velbinsize=8;%velocity bin size (n x n square)
gradbinsize=16;%gradient bin size (n x n square)
edgecutoff=1;%bins on the edge cutoff

%PTV parameters
particlediameter=6;%maximum pixel diameter of tracer particle
minintensity=40;%minimum tracer particle intensity to be considered a feature, change this so between 300 and 500 features are found for a well conditioned system
searchrad=8;%pixel radius to search for particle between frames
frameappear=5;%frames a particle must appear in to be considered a feature
framedisappear=5;%frames a particle can disappear for before forgotten about



%---------------------------------------------------
tic
disp(strcat('Reading Images...'));
%read images from file
r1=[];
for i=1:numbofim
    im=imread(strcat(imageloc,'\',imagename,num2str(i-1,'%04.0f'),'.tif'));
    r=Nfeature3(im,1,particlediameter,2,minintensity); %(file,lengthscale of noise,diameter of feature,2D,minimum intensity)
%     r(:,2)=480-r(:,2);
    temp=r;
    temp(:,6)=i;
    r1=[r1;temp];
end
toc
disp(strcat('Tracking Particles...'));

%ptv calculation
s1=r1(:,1:2);%organize position data for trackmem
s1(:,3)=r1(:,6);
[lub1]=trackmem(s1,searchrad,2,frameappear,framedisappear);%(xyzs,set blocksize,'2'D,frames feature appears,frames between reappear)
toc

lub=[];
lub=[lub1];

disp(strcat('Extracting Velocity Field...'));

%extract velocity information from particle tracking
velocity=[];
df=[];

for i=1:numbofim
    start1=lub(find(lub(:,3)==i),:);
    finish1=lub(find(lub(:,3)==i+trajectorygap),:);
    [id ia ib]=intersect(start1(:,4),finish1(:,4),'rows');
    pt11=start1(ia,:);
    pt12=finish1(ib,:);
    a=size(pt11);
    for j=1:a(1)
      df(j,1)=pt11(j,1);
      df(j,2)=pt11(j,2);
      df(j,3)=(pt12(j,1)-pt11(j,1))/trajectorygap;
      df(j,4)=(pt12(j,2)-pt11(j,2))/trajectorygap;
    end   
    velocity=[velocity;df];
end
toc
disp(strcat('Binning Velocity Field...'));

%bin velocity
[binvelocity,matbinvelu,matvinvelv] = fbinvel(velocity,imagexres,imageyres,velbinsize,edgecutoff);%change velocity matrix into binned matrixes

%plot velocity magnitude
s=velbinsize;
x=[];
y=[];
c=[];
x=binvelocity(:,1)*lengthperpix-xcenter*lengthperpix;
y=binvelocity(:,2)*lengthperpix-ycenter*lengthperpix;
c=sqrt((binvelocity(:,3).^2+binvelocity(:,4).^2))*lengthperpix*framerate;

figure
colormap jet;
scatter(x,y,s,c,'filled','s'),colorbar;axis equal;axis([(0-xcenter)*lengthperpix (imagexres-xcenter)*lengthperpix (0-ycenter)*lengthperpix (imageyres-ycenter)*lengthperpix]);
xlabel(lengthunit)
ylabel(lengthunit)
title(colorbar,strcat('Velocity Magnitude (',lengthunit,'/',timeunit,')'))


%ends velocity field calculation
 
%begin gradient calculation 
disp(strcat('Calculating Gradients in x Direction...')); 
dpux=(imagexres*imageyres/(yres*yres),3);
dpvx=zeros(imagexres*imageyres/(yres*yres),3);

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
        pux=polyfit(tempdx(:,1),tempdx(:,3)*framerate,polyfitdeg); %fit polynomial to u
        pvx=polyfit(tempdx(:,1),tempdx(:,4)*framerate,polyfitdeg); %fit polynomial to v
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
dpuy=zeros(imagexres*imageyres/(xres*xres),3);
dpvy=zeros(imagexres*imageyres/(xres*xres),3);
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
            dpuy(o,1)=avgx;
            dpuy(o,2)=y;
            dpuy(o,3)=polyval(puy,y);
            o=o+1;
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
bindpux=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
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
bindpvx=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
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
bindpuy=zeros(imagexres*imageyres/(xres*xres*gradbinsize*gradbinsize),3);
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
bindpvy=zeros(imagexres*imageyres/(xres*xres*gradbinsize*gradbinsize),3);
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
scatter((FTP(:,1)*gradbinsize-xcenter)*lengthperpix,(FTP(:,2)*gradbinsize-ycenter)*lengthperpix,0.5*gradbinsize^2,FTP(:,3),'filled','s')
caxis([-1,1])
colorbar;axis equal;axis([(roixmin-xcenter)*lengthperpix (roixmax-xcenter)*lengthperpix (roiymin-ycenter)*lengthperpix (roiymax-ycenter)*lengthperpix]);
xlabel(lengthunit)
ylabel(lengthunit)
title(colorbar,strcat('FTP'))

figure
colormap jet;
scatter((MVGT(:,1)*gradbinsize-xcenter)*lengthperpix,(MVGT(:,2)*gradbinsize-ycenter)*lengthperpix,0.5*gradbinsize^2,MVGT(:,3),'filled','s')
caxis([0,maxMVGT])
colorbar;axis equal;axis([(roixmin-xcenter)*lengthperpix (roixmax-xcenter)*lengthperpix (roiymin-ycenter)*lengthperpix (roiymax-ycenter)*lengthperpix]);
xlabel(lengthunit)
ylabel(lengthunit)
title(colorbar,strcat('MVGT (1/',timeunit,')'))
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

save(strcat(imageloc,'\',imagename,'matlabdata.mat'))%save generated matlab data