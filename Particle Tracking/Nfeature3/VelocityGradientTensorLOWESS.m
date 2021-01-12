%Analyze set of processed images for velocity, FTP, and MVGT fields
%functions required:
%Nfeature3.m (and functions that this one uses)
%trackmem.m (and functions that this one uses)

%inputs---------------------------------------------

%image parameters
lengthperpix=5/1058;%distance per pixel
lengthunit='mm'; %unit of distance
framerate=30;%frames per time
timeunit='sec'; %unit of time
imagexres=2048;%pixel resolution of images
imageyres=1088;
numbofim=200;%number of images
imageloc='C:\Users\Patrick\16-10-12\Glycerol300ppm 30fps mFoRMv2 Q2=1.0mlmin\FRR=0.8';%image files location
imagename='1.0mlminfrr0.830fps(2)';%image name without '0000','0001',etc or '.tif'

%graphing parameters
xcenter=954;%pixel center of geometry
ycenter=520;
maxMVGT=1;%colorbar maximum for MVGT (1/timeunit)

%velocity and gradient quality parameters
trajectorygap=2; %images between velocity calculation
p=30;
roixmin=400;
roixmax=1500;
roiymin=15;
roiymax=1065;

%binning parameters
velbinsize=8;%velocity bin size (n x n square)
gradbinsize=8;%gradient bin size (n x n square)
edgecutoff=1;%bins on the edge cutoff

%PTV parameters
particlediameter=6;%maximum pixel diameter of tracer particle
minintensity=40;%minimum tracer particle intensity to be considered a feature
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
        if ((pt11(j,1)+pt12(j,1))/2>roixmin) && ((pt11(j,1)+pt12(j,1))/2<roixmax) && ((pt11(j,2)+pt12(j,2))/2>roiymin) && ((pt11(j,2)+pt12(j,2))/2<roiymax) %only evaluate velocity in the region of interest
            df(j,1)=(pt11(j,1)+pt12(j,1))/2;
            df(j,2)=(pt11(j,2)+pt12(j,2))/2;
            df(j,3)=(pt12(j,1)-pt11(j,1))/trajectorygap;
            df(j,4)=(pt12(j,2)-pt11(j,2))/trajectorygap;
        end
    end
    df(any(df==0,2),:)=[];
    velocity=[velocity;df];
end
toc
disp(strcat('Binning Velocity Field...'));

%bin velocity
[binvelocity,matbinvelu,matvinvelv] = fbinvel(velocity,imagexres,imageyres,velbinsize,edgecutoff);%change velocity matrix into binned matrixes

%plot velocity magnitude
s=velbinsize*4;
x=[];
y=[];
c=[];
x=binvelocity(:,1)*lengthperpix-xcenter*lengthperpix;
y=binvelocity(:,2)*lengthperpix-ycenter*lengthperpix;
c=sqrt((binvelocity(:,3).^2+binvelocity(:,4).^2))*lengthperpix*framerate;

figure
colormap jet;
scatter(x,y,s,c,'filled','s'),colorbar;axis equal;axis([(roixmin-xcenter)*lengthperpix (roixmax-xcenter)*lengthperpix (roiymin-ycenter)*lengthperpix (roiymax-ycenter)*lengthperpix]);
xlabel(lengthunit)
ylabel(lengthunit)
title(colorbar,strcat('Velocity Magnitude (',lengthunit,'/',timeunit,')'))


%ends velocity field calculation
 
%begin gradient calculation 
disp(strcat('Calculating Gradients...')); 

a=size(velocity);
dpux=zeros(a(1),4);
dpvx=zeros(a(1),4);
dpuy=zeros(a(1),4);
dpvy=zeros(a(1),4);
fitstat=zeros(a(1),2);
fitmat=zeros(a(1),5);
tic
for i=1:10:a(1)
    g=1;
    for m=1:a(1)
        r=sqrt((velocity(i,1)-velocity(m,1))^2+(velocity(i,2)-velocity(m,2))^2);
        if(r<p)        
            fitmat(g,1)=velocity(m,1);
            fitmat(g,2)=velocity(m,2);
            fitmat(g,3)=velocity(m,3);
            fitmat(g,4)=velocity(m,4);
            fitmat(g,5)=(1-(r/p)^3)^3;
            g=g+1;
        end
    end
    fitmat(any(fitmat==0,2),:)=[];
    [fu,gofu]=fit([(fitmat(:,1)-velocity(i,1)) (fitmat(:,2)-velocity(i,2))],fitmat(:,3),'poly22','Weights',fitmat(:,5));
    [fv,gofv]=fit([(fitmat(:,1)-velocity(i,1)) (fitmat(:,2)-velocity(i,2))],fitmat(:,4),'poly22','Weights',fitmat(:,5));
    ucoeff=coeffvalues(fu);
    vcoeff=coeffvalues(fv);
    dpux(i,1)=velocity(i,1);
    dpuy(i,1)=velocity(i,1);
    dpvx(i,1)=velocity(i,1);
    dpvy(i,1)=velocity(i,1);
    dpux(i,2)=velocity(i,2);
    dpuy(i,2)=velocity(i,2);
    dpvx(i,2)=velocity(i,2);
    dpvy(i,2)=velocity(i,2);
    dpux(i,3)=ucoeff(2)*framerate;
    dpuy(i,3)=ucoeff(3)*framerate;
    dpvx(i,3)=vcoeff(2)*framerate;
    dpvy(i,3)=vcoeff(3)*framerate;
    dpux(i,4)=(1/gofu.rmse)^2;
    dpuy(i,4)=(1/gofu.rmse)^2;
    dpvx(i,4)=(1/gofv.rmse)^2;
    dpvy(i,4)=(1/gofv.rmse)^2;
    if (mod(i-1,1000)==0)
        complete=round(100*i/a(1),1)
    end
end
toc


%remove zeros from matrix
dpux(any(dpux==0,2),:)=[];
dpvx(any(dpvx==0,2),:)=[];
dpuy(any(dpuy==0,2),:)=[];
dpvy(any(dpvy==0,2),:)=[];

disp(strcat('Binning du/dx...'));

%bin du/dx
a=size(dpux);
tic
matbindpux=zeros(imagexres/(gradbinsize),imageyres/(gradbinsize));
bindpux=zeros(imagexres*imageyres/(gradbinsize*gradbinsize),3);
counter=zeros(imagexres/(gradbinsize),imageyres/(gradbinsize));
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
matbindpvx=zeros(imagexres/(gradbinsize),imageyres/(gradbinsize));
bindpvx=zeros(imagexres*imageyres/(gradbinsize*gradbinsize),3);
counter=zeros(imagexres/(gradbinsize),imageyres/(gradbinsize));
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
matbindpuy=zeros(imagexres/(gradbinsize),imageyres/(gradbinsize));
bindpuy=zeros(imagexres*imageyres/(gradbinsize*gradbinsize),3);
counter=zeros(imagexres/(gradbinsize),imageyres/(gradbinsize));
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
matbindpvy=zeros(imagexres/(gradbinsize),imageyres/(gradbinsize));
bindpvy=zeros(imagexres*imageyres/(gradbinsize*gradbinsize),3);
counter=zeros(imagexres/(gradbinsize),imageyres/(gradbinsize));
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
continuity=zeros(imagexres/(gradbinsize),imageyres/(gradbinsize));
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
scatter((FTP(:,1)*gradbinsize-xcenter)*lengthperpix,(FTP(:,2)*gradbinsize-ycenter)*lengthperpix,gradbinsize*2,FTP(:,3),'filled','s')
caxis([-1,1])
colorbar;axis equal;axis([(roixmin-xcenter)*lengthperpix (roixmax-xcenter)*lengthperpix (roiymin-ycenter)*lengthperpix (roiymax-ycenter)*lengthperpix]);
xlabel(lengthunit)
ylabel(lengthunit)
title(colorbar,strcat('FTP'))

figure
colormap jet;
scatter((MVGT(:,1)*gradbinsize-xcenter)*lengthperpix,(MVGT(:,2)*gradbinsize-ycenter)*lengthperpix,gradbinsize*2,MVGT(:,3),'filled','s')
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