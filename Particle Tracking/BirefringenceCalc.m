clear all
X=zeros(1,4);

%Load images
imageloc0='C:\Users\Patrick\Microfourrollmill CNC\Flow Type\300 degrees';%image files location
imagename0='lambda1.0V1.0';%image name without '0000','0001',etc or '.tif'
%firstim0 = 18; %first frame number
%lastim0 = 800; %last frame number

imageloc45='C:\Users\Patrick\Microfourrollmill CNC\Flow Type\345 degrees';%image files location
imagename45='lambda1.0V1.0';%image name without '0000','0001',etc or '.tif'
%firstim45 = 10; %first frame number
%lastim45 = 800; %last frame number

imageloc_0='C:\Users\Patrick\Microfourrollmill CNC\Flow Type\300 degrees';%image files location
imagename_0='maximum intensity';%image name without '0000','0001',etc or '.tif'

imageloc_bkgd='C:\Users\Patrick\Microfourrollmill CNC\Flow Type\300 degrees';%image files location
imagename_bkgd='background';%image name without '0000','0001',etc or '.tif'

%image segment
xmin = 930;
xmax = 1080;
ymin = 350;
ymax = 500;

numbofim = 1;%min(lastim0-firstim0, lastim45-firstim45); %number of images

delta_all = zeros(1,numbofim);
theta_all = zeros(1,numbofim);

tic
disp(strcat('Reading Images...'));
%read images from file
for i= 1:numbofim
    im0_uint8=imread(strcat(imageloc0,'\',imagename0,'.tif'));%,num2str(firstim0+i-1,'%04.0f'),'.tif'));
    im0 = double(im0_uint8);
    im45_uint8=imread(strcat(imageloc45,'\',imagename45,'.tif'));%num2str(firstim45+i-1,'%04.0f'),'.tif'));
    im45 = double(im45_uint8);
    im_0_uint8=imread(strcat(imageloc_0,'\',imagename_0,'.tif'));%num2str(310,'%04.0f'),'.tif'));
    im_0 = double(im_0_uint8);
    im_bkgd_uint8=imread(strcat(imageloc_bkgd,'\',imagename_bkgd,'.tif'));%,num2str(310,'%04.0f'),'.tif'));
    im_bkgd = double(im_bkgd_uint8);
    
    i0 = 0.5 * (im0-im_bkgd) ./ im_0;   %Factor 2 without any cross polarizers, 1 with one and 0.5 with 2 cross polarizers
    i45 = 0.5 * (im45-im_bkgd) ./ im_0; %Factor 2 without any cross polarizers, 1 with one and 0.5 with 2 cross polarizers
    
    if i0 == 0;
        i0 = 10^(-9);
    end
    if i45 == 0;
        i45 = 10^(-9);
    end
    
    delta = asin(sqrt(i0+i45));
    %delta_uint8 = uint8(delta/pi*360);
    meandelta = real(mean(mean(delta(ymin:ymax, xmin:xmax)),2)); %average of a segment of the image. First y value then x value
    delta_all(i) = meandelta;
    av_delta = nanmean(delta_all);
    std_delta = std(delta_all);
    
    theta = 0.5 * atan(sqrt(i0./i45));
    %theta_uint8 = uint8(theta/pi*360);
    meantheta = real(mean(mean(theta(ymin:ymax, xmin:xmax)),2));
    theta_all(i) = meantheta;
    av_theta = nanmean(theta_all);
    std_theta = nanstd(theta_all);
    X=[X;av_delta   std_delta  av_theta  std_theta];
end
toc
% 
%plot(theta_all)
% saveas(figure(1), 'C:\Users\Patrick\18-8-14 (Wormlike micelles)\Data\Q2=11mlminQ1=1.1mlmin10fps_theta.fig');
% close
%plot(delta_all)
% saveas(figure(1), 'C:\Users\Patrick\18-8-14 (Wormlike micelles)\Data\Q2=11mlminQ1=1.1mlmin10fpsdelta.fig');
% save('C:\Users\Patrick\18-7-23\5% CNC with 10%Glycerol\Data\Q2=11mlminQ1=1.1mlmin10fps.mat', 'delta_all');


% xfit = (1:59)';
% yfit = (delta_all(360:end))';
% 
% f = @(b,x) b(1).*exp(-x./b(2)) + b(3);
% fit = fminsearch(@(b) norm(yfit - f(b,xfit)), [;5;0.1]);
%  
% % delta_max = mean(delta_all(40:90));
% delta_max = max(delta_all);
% 
% save('C:\Users\Patrick\18-7-23\5% CNC with 10%Glycerol\Data\Q2=0.07mlminQ1=0.07mlminfps15_rel_30s.mat', 'fit', 'delta_max', 'delta_all');
