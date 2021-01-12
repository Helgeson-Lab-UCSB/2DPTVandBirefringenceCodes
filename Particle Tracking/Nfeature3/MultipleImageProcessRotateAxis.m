r1=[]
for i=1:59
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Calibration with Karo Light Corn Syrup\2D Extensional Flow Stepsize 100000\Measurement1\Measurement10',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,40);
    temp=r;
    temp(:,6)=i;
    r1=[r1;temp];
end

r2=[]
for i=1:59
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Calibration with Karo Light Corn Syrup\2D Extensional Flow Stepsize 100000\Measurement2\Measurement20',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,40);
    temp=r;
    temp(:,6)=i;
    r2=[r2;temp];
end

r3=[]
for i=1:59
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Calibration with Karo Light Corn Syrup\2D Extensional Flow Stepsize 100000\Measurement3\Measurement30',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,40);
    temp=r;
    temp(:,6)=i;
    r3=[r3;temp];
end

r4=[]
for i=1:59
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Calibration with Karo Light Corn Syrup\2D Extensional Flow Stepsize 100000\Measurement4\Measurement40',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,40);
    temp=r;
    temp(:,6)=i;
    r4=[r4;temp];
end



s1=r1(:,1:2);
s1(:,3)=r1(:,6);

s2=r2(:,1:2);
s2(:,3)=r2(:,6);

s3=r3(:,1:2);
s3(:,3)=r3(:,6);

s4=r4(:,1:2);
s4(:,3)=r4(:,6);


[lub1]=trackmem(s1,5,2,5,5);
[lub2]=trackmem(s2,5,2,5,5);
[lub3]=trackmem(s3,5,2,5,5);
[lub4]=trackmem(s4,5,2,5,5);

length1=max(lub1(:,4));
length2=max(lub2(:,4));
length3=max(lub3(:,4));
length4=max(lub4(:,4));

lub2(:,4)=lub2(:,4)+length1;
lub3(:,4)=lub3(:,4)+length1+length2;
lub4(:,4)=lub4(:,4)+length1+length2+length3;
    
lub=[];
lub=[lub1;lub2;lub3;lub4];

% totallength=length1+length2+length3+length4;
% for k=1:totallength
%     pt=lub(find(lub(:,4)==k),:);
%     pt=sortrows(pt,3);
%     hold on;
%     if k<=length1
%         plot(pt(:,1),pt(:,2),'b')
%     elseif k<=length1+length2
%            plot(pt(:,1),pt(:,2),'r')  
%     elseif k<=length1+length2+length3
%            plot(pt(:,1),pt(:,2),'g')
%     else plot(pt(:,1),pt(:,2),'k')
%     end
% end
% box on;

start1=lub(find(lub(:,3)==1),:);
finish1=lub(find(lub(:,3)==15),:);
[id ia ib]=intersect(start1(:,4),finish1(:,4),'rows');
pt11=start1(ia,:);
pt12=finish1(ib,:);
a=size(pt11);

% start2=lub(find(lub(:,3)==16),:);
% finish2=lub(find(lub(:,3)==30),:);
% [id ia ib]=intersect(start2(:,4),finish2(:,4),'rows');
% pt21=start2(ia,:);
% pt22=finish2(ib,:);
% b=size(pt21);
% 
% start3=lub(find(lub(:,3)==31),:);
% finish3=lub(find(lub(:,3)==45),:);
% [id ia ib]=intersect(start3(:,4),finish3(:,4),'rows');
% pt31=start3(ia,:);
% pt32=finish3(ib,:);
% c=size(pt31);
% 
%  looplength=a(1)+b(1)+c(1);
%  figure
%  for m=1:looplength
%      if m<=a(1)
%          vectarrow(pt11(m,1:2),pt12(m,1:2));
%      elseif m<=a(1)+b(1)
%          vectarrow(pt21(m-a(1),1:2),pt22(m-a(1),1:2));
%      else
%          vectarrow(pt31(m-a(1)-b(1),1:2),pt32(m-a(1)-b(1),1:2));
%      end
%      hold on;
%  end
% % 

xdf=[];
ydf=[];
xdf(:,1)=pt11(:,1);
xdf(:,2)=pt11(:,2);
xdf(:,3)=pt12(:,1)-pt11(:,1);
ydf(:,1)=pt11(:,1);
ydf(:,2)=pt11(:,2);
ydf(:,3)=pt12(:,2)-pt11(:,2);

% for t=1:a(1)
%     plot(xdf(t,1),xdf(t,2),'k')
%     plot(ydf(t,1),ydf(t,2),'r')
%     hold on;
% end
 
xcenterline=[];
ycenterline=[];
v=1;
w=1;
eta=5;
 for u=1:a(1)
     if abs(xdf(u,3))<=eta
         xcenterline(v,1)=xdf(u,1);
         xcenterline(v,2)=xdf(u,2);
         v=v+1;
     end
     if abs(ydf(u,3))<=eta
         ycenterline(w,1)=ydf(u,1);
         ycenterline(w,2)=ydf(u,2);
         w=w+1;
     end
 end
 
b=size(xcenterline);
c=size(ycenterline);
 
 figure
  for t=1:b(1)
     plot(xcenterline(t,1),xcenterline(t,2),'k')
     hold on;
  end
  
  figure
  for p=1:c(1)
      plot(ycenterline(p,1),ycenterline(p,2),'r')
      hold on;
  end
%  xdf(:,1)=abs(xdf(:,1)-xcenter);
%  xdf(:,2)=abs(xdf(:,2));
%  ydf(:,1)=abs(ydf(:,1)-ycenter);
%  ydf(:,2)=abs(ydf(:,2));