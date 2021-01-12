%Flow type parameter (FTP) and magnitude of velocity gradient tensor (MVGT) calculation
MagE=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
MagO=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
FTP=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
MVGT=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
G=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
stretch=zeros(imagexres*imageyres/(yres*yres*gradbinsize*gradbinsize),3);
continuity=zeros(imagexres/(xres*gradbinsize),imageyres/(xres*gradbinsize));
lambda=0.0569;
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
            G(v,1)=i;
            G(v,2)=k;
            G(v,3)=2.*MVGT(v,3)./sqrt(1+(FTP(v,3)).^2);
            stretch(v,1)=i;
            stretch(v,2)=k;
            stretch(v,3)=G(v,3).*sqrt(FTP(v,3)).*lambda;
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
G(any(G==0,2),:)=[];
stretch(any(stretch==0,2),:)=[];

%graph FTP, MVGT, and extensional rate
figure(1)
colormap jet;
scatter((FTP(:,1)*gradbinsize-xcenter)*lengthperpix,(FTP(:,2)*gradbinsize-ycenter)*lengthperpix,gradbinsize*4,FTP(:,3),'filled','s')
caxis([-1,1])
colorbar;axis equal;axis([(0-xcenter)*lengthperpix (imagexres-xcenter)*lengthperpix (0-ycenter)*lengthperpix (imageyres-ycenter)*lengthperpix]);
xlabel(lengthunit)
ylabel(lengthunit)
title(colorbar,strcat('FTP'))
xlim([-2,2])
ylim([-2,2])
axis square

figure(2)
colormap jet;
scatter((stretch(:,1)*gradbinsize-xcenter)*lengthperpix,(stretch(:,2)*gradbinsize-ycenter)*lengthperpix,gradbinsize*4,stretch(:,3),'filled','s')
colorbar;axis equal;axis([(0-xcenter)*lengthperpix (imagexres-xcenter)*lengthperpix (0-ycenter)*lengthperpix (imageyres-ycenter)*lengthperpix]);
caxis([0,1])
xlabel(lengthunit)
ylabel(lengthunit)
title(colorbar,strcat('Wi sqrt(\Lambda)'))
xlim([-2,2])
ylim([-2,2])
axis square

%Calculate Average Extensional Rate within some radius of the center
beamsize=0.5/2;%Radius of area of interest in units of length unit

counter=0;
AVGFTP=0;
AVGstretch=0;

%Calculate Average FTP within some radius of the center
counter=0;
AVGFTP=0;
v=length(FTP(:,1));
for i=1:v-1
    if sqrt((FTP(i,1)*gradbinsize-xcenter)*lengthperpix*(FTP(i,1)*gradbinsize-xcenter)*lengthperpix+(FTP(i,2)*gradbinsize-ycenter)*lengthperpix*(FTP(i,2)*gradbinsize-ycenter)*lengthperpix)<beamsize
        AVGFTP=FTP(i,3)+AVGFTP;
        counter=counter+1;
    end    
end
AVGFTP=AVGFTP/counter;

%Calculate Variance in FTP in that radius
counter=0;
STDFTP=0;
for i=1:v-1
    if sqrt((FTP(i,1)*gradbinsize-xcenter)*lengthperpix*(FTP(i,1)*gradbinsize-xcenter)*lengthperpix+(FTP(i,2)*gradbinsize-ycenter)*lengthperpix*(FTP(i,2)*gradbinsize-ycenter)*lengthperpix)<beamsize
        STDFTP=(FTP(i,3)-AVGFTP).^2+STDFTP;
        counter=counter+1;
    end    
end
STDFTP=sqrt(STDFTP/counter);

%Calculate Average stretch within some radius of the center
counter=0;
AVGstretch=0;
v=length(G(:,1));
for i=1:v-1
    if sqrt((stretch(i,1)*gradbinsize-xcenter)*lengthperpix*(stretch(i,1)*gradbinsize-xcenter)*lengthperpix+(stretch(i,2)*gradbinsize-ycenter)*lengthperpix*(stretch(i,2)*gradbinsize-ycenter)*lengthperpix)<beamsize
        AVGstretch=stretch(i,3)+AVGstretch;
        counter=counter+1;
    end    
end
AVGstretch=AVGstretch/counter;

%Calculate Variance in G in that radius
counter=0;
STDstretch=0;
for i=1:v-1
    if sqrt((stretch(i,1)*gradbinsize-xcenter)*lengthperpix*(stretch(i,1)*gradbinsize-xcenter)*lengthperpix+(stretch(i,2)*gradbinsize-ycenter)*lengthperpix*(stretch(i,2)*gradbinsize-ycenter)*lengthperpix)<beamsize
        STDstretch=(stretch(i,3)-AVGstretch).^2+STDstretch;
        counter=counter+1;
    end    
end
STDstretch=sqrt(STDstretch/counter);

X=[real(AVGstretch) real(STDstretch)  real(AVGFTP)  real(STDFTP)];