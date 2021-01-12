%Calculate Average Extensional Rate within some radius of the center
beamsize=0.5/2;%Radius of area of interest in units of length unit

counter=0;
AVGFTP=0;
AVGEXT=0;
% v=length(bindpux(:,1));
% for i=1:v-1
%     if sqrt((bindpux(i,1)-xcenter)*lengthperpix*(bindpux(i,1)-xcenter)*lengthperpix+(bindpux(i,2)-ycenter)*lengthperpix*(bindpux(i,2)-ycenter)*lengthperpix)<beamsize
%         AVGEXT=bindpux(i,3)+AVGEXT;
%         counter=counter+1;
%     end    
% end
% 
% AVGEXT=AVGEXT/counter
% 
% %Calculate Variance in Average Extensional Rate
% counter=0;
% VAREXT=0;
% for i=1:v-1
%     if sqrt((bindpux(i,1)-xcenter)*lengthperpix*(bindpux(i,1)-xcenter)*lengthperpix+(bindpux(i,2)-ycenter)*lengthperpix*(bindpux(i,2)-ycenter)*lengthperpix)<beamsize
%         VAREXT=(bindpux(i,3)-AVGEXT).^2+VAREXT;
%         counter=counter+1;
%     end    
% end
% VAREXT=VAREXT/counter

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

%Calculate Average MVGT within some radius of the center
counter=0;
AVGMVGT=0;
v=length(MVGT(:,1));
for i=1:v-1
    if sqrt((MVGT(i,1)*gradbinsize-xcenter)*lengthperpix*(MVGT(i,1)*gradbinsize-xcenter)*lengthperpix+(MVGT(i,2)*gradbinsize-ycenter)*lengthperpix*(MVGT(i,2)*gradbinsize-ycenter)*lengthperpix)<beamsize
        AVGMVGT=MVGT(i,3)+AVGMVGT;
        counter=counter+1;
    end    
end
AVGMVGT=AVGMVGT/counter;

%Calculate Variance in MVGT in that radius
counter=0;
STDMVGT=0;
for i=1:v-1
    if sqrt((MVGT(i,1)*gradbinsize-xcenter)*lengthperpix*(MVGT(i,1)*gradbinsize-xcenter)*lengthperpix+(MVGT(i,2)*gradbinsize-ycenter)*lengthperpix*(MVGT(i,2)*gradbinsize-ycenter)*lengthperpix)<beamsize
        STDMVGT=(MVGT(i,3)-AVGMVGT).^2+STDMVGT;
        counter=counter+1;
    end    
end
STDMVGT=sqrt(STDMVGT/counter);

%Calculate Average G within some radius of the center
counter=0;
AVGG=0;
v=length(G(:,1));
for i=1:v-1
    if sqrt((G(i,1)*gradbinsize-xcenter)*lengthperpix*(G(i,1)*gradbinsize-xcenter)*lengthperpix+(G(i,2)*gradbinsize-ycenter)*lengthperpix*(G(i,2)*gradbinsize-ycenter)*lengthperpix)<beamsize
        AVGG=G(i,3)+AVGG;
        counter=counter+1;
    end    
end
AVGG=AVGG/counter;

%Calculate Variance in G in that radius
counter=0;
STDG=0;
for i=1:v-1
    if sqrt((G(i,1)*gradbinsize-xcenter)*lengthperpix*(G(i,1)*gradbinsize-xcenter)*lengthperpix+(G(i,2)*gradbinsize-ycenter)*lengthperpix*(G(i,2)*gradbinsize-ycenter)*lengthperpix)<beamsize
        STDG=(G(i,3)-AVGG).^2+STDG;
        counter=counter+1;
    end    
end
STDG=sqrt(STDG/counter);

%Calculate Average E within some radius of the center
counter=0;
AVGE=0;
v=length(MagE(:,1));
for i=1:v-1
    if sqrt((MagE(i,1)*gradbinsize-xcenter)*lengthperpix*(MagE(i,1)*gradbinsize-xcenter)*lengthperpix+(MagE(i,2)*gradbinsize-ycenter)*lengthperpix*(MagE(i,2)*gradbinsize-ycenter)*lengthperpix)<beamsize
        AVGE=MagE(i,3)+AVGE;
        counter=counter+1;
    end    
end
AVGE=AVGE/counter;

X=[X;AVGG   STDG  AVGFTP  STDFTP AVGMVGT STDMVGT AVGE];