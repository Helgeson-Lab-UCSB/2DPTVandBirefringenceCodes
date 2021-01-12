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