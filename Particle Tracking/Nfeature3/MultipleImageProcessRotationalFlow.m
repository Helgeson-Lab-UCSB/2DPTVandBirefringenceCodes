r1=[]
for i=1:59
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Measurement1\Measurement10',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r1=[r1;temp];
end

r2=[]
for i=1:59
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Measurement2\Measurement20',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp(:,6)=i;
    r2=[r2;temp];
end

r3=[]
for i=1:59
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Measurement3\Measurement30',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r3=[r3;temp];
end

r4=[]
for i=1:59
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Measurement4\Measurement40',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
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

% totallength=size(lub);
% for k=1:totallength
%     pt=lub(find(lub(:,4)==k),:);
%     pt=sortrows(pt,3);
%     hold on;
%     plot(pt(:,1),pt(:,2),'b')
% end
% box on;

velocity=[];
df=[];
for i=1:59
    start1=lub(find(lub(:,3)==i),:);
    finish1=lub(find(lub(:,3)==i+1),:);
    [id ia ib]=intersect(start1(:,4),finish1(:,4),'rows');
    pt11=start1(ia,:);
    pt12=finish1(ib,:);
    a=size(pt11);
    for j=1:a(1)
        df(j,1)=pt11(j,1);
        df(j,2)=pt11(j,2);
        df(j,3)=pt12(j,1)-pt11(j,1);
        df(j,4)=pt12(j,2)-pt11(j,2);
    end
    velocity=[velocity;df];
end

a=size(velocity);
%  figure
%   for t=1:a(1)
%      plot(velocity(t,1),velocity(t,3),'k')
%      plot(velocity(t,2),velocity(t,4),'r')
%      hold on;
%  end

binvelocity=[];
t=1;
for i=0:63
    for k=0:47
        bin=[];
        for j=1:a(1)
            if and(velocity(j,1)>=i*10,velocity(j,1)<=(i+1)*10)&& and(velocity(j,2)>=k*10,velocity(j,2)<=(k+1)*10)
                    bin=[bin;velocity(j,:)];
            end
        end
        y=size(bin);
        if y(1)>0
            binvelocity(t,1)=(2*i+1)/2;
            binvelocity(t,2)=(2*k+1)/2;
            binvelocity(t,3)=mean(bin(:,3));
            binvelocity(t,4)=mean(bin(:,4));
            t=t+1;
        end
    end
end
% 
b=size(binvelocity);
% xnormalvelocity=[];
% for i=0:63
%     temp=[];
%     for j=1:b(1)
%         if and(binvelocity(j,1)<i+1, binvelocity(j,1)>i)
%             temp=[temp;binvelocity(j,:)];
%         end
%     end
%     xnormalvelocity(i+1,1)=temp(1,1);
%     xnormalvelocity(i+1,2)=mean(temp(:,3));
%     xnormalvelocity(i+1,3)=std(temp(:,3));
%     xnormalvelocity(i+1,4)=abs(std(temp(:,3))/mean(temp(:,3)));
% end
% 
% ynormalvelocity=[];
% for i=0:47
%     temp=[];
%     for j=1:b(1)
%         if and(binvelocity(j,2)<i+1, binvelocity(j,2)>i)
%             temp=[temp;binvelocity(j,:)];
%         end
%     end
%     ynormalvelocity(i+1,1)=temp(1,2);
%     ynormalvelocity(i+1,2)=mean(temp(:,4));
%     ynormalvelocity(i+1,3)=std(temp(:,4));
%     ynormalvelocity(i+1,4)=abs(std(temp(:,3))/mean(temp(:,3)));
% end



% 
% 
%  looplength=size(velocity);
%  figure
%  for m=1:looplength(1)
%      vectarrow(velocity(m,1:2),velocity(m,3:4));
%      hold on;
%  end
% % % 
% 
% 
% 
% for t=1:a(1)
%     plot(xdf(t,1),xdf(t,2),'k')
%     plot(ydf(t,1),ydf(t,2),'r')
%     hold on;
% end
%  a=size(velocity);
 xcenter=binvelocity(1,1);
 ycenter=binvelocity(1,2);
 xinitialerror=abs(binvelocity(1,4)-0);
 yinitialerror=abs(binvelocity(1,3)-0);
 for u=2:b(1)
     xerror=abs(binvelocity(u,4)-0);
     yerror=abs(binvelocity(u,3)-0);
     if (xerror<xinitialerror)&& and (binvelocity(u,1)>=20,binvelocity(u,1)<=40)
         xcenter=binvelocity(u,1);
         xinitialerror=xerror;
     else 
         xcenter=xcenter;   
     end
     if (yerror<yinitialerror) && and (binvelocity(u,2)>=15,binvelocity(u,2)<=33)
         ycenter=binvelocity(u,2);
         yinitialerror=yerror;
     else 
         ycenter=ycenter;   
     end
 end
%

rotationvelocity=[];
for j=1:b(1)
    rotationvelocity(j,1)=sqrt((binvelocity(j,1)-xcenter)^2+(binvelocity(j,2)-ycenter)^2);
    rotationvelocity(j,2)=sqrt((binvelocity(j,3))^2+(binvelocity(j,4))^2);
end
%  velocity(:,1)=abs(velocity(:,1)-xcenter);
%  velocity(:,2)=abs(velocity(:,2)-ycenter);
%  velocity(:,3)=abs(velocity(:,3));
%  velocity(:,4)=abs(velocity(:,4));
%  
newbinvelocity=[]; 
c=size(rotationvelocity);
t=1;
for i=0:60
    bin=[];
    for j=1:c(1)
        if and(rotationvelocity(j,1)>=i,rotationvelocity(j,1)<=(i+1))
                    bin=[bin;rotationvelocity(j,:)];
        end
    end
    y=size(bin);
    if y(1)>0
       newbinvelocity(t,1)=(2*i+1)/2;
       newbinvelocity(t,2)=mean(bin(:,2));
       t=t+1;
    end
 end