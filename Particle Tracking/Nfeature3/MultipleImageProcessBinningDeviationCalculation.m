r1=[]
for i=1:59
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Measurement1\Measurement10',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,80);
    temp=r;
    temp(:,6)=i;
    r1=[r1;temp];
end

r2=[]
for i=1:59
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Measurement2\Measurement20',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,80);
    temp=r;
    temp(:,6)=i;
    r2=[r2;temp];
end

r3=[]
for i=1:59
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Measurement3\Measurement30',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,80);
    temp=r;
    temp(:,6)=i;
    r3=[r3;temp];
end

r4=[]
for i=1:59
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Measurement4\Measurement40',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,80);
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
shearrate=1;
xshearrateintercept=0.2;
yshearrateintercept=0.5;

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
            binvelocity(t,1)=(2*i+1)/2*10;
            binvelocity(t,2)=(2*k+1)/2*10;
            binvelocity(t,3)=mean(bin(:,3))*15;
            binvelocity(t,4)=mean(bin(:,4))*15;
            nominalxvelocity=-shearrate*binvelocity(t,1)+xshearrateintercept;
            nominalyvelocity=shearrate*binvelocity(t,2)+yshearrateintercept;
            binvelocity(t,5)=binvelocity(t,3)-nominalxvelocity;
            binvelocity(t,6)=binvelocity(t,4)-nominalyvelocity;

%             if binvelocity(t,6)>=100
%                 binvelocity(t,6)=100;
%             end
            t=t+1;
        end
    end
end

b=size(binvelocity);
xnormalvelocity=[];
for i=0:63
    temp=[];
    for j=1:b(1)
        if and(binvelocity(j,1)<(i+1)*10, binvelocity(j,1)>i*10)
            temp=[temp;binvelocity(j,:)];
        end
    end
    si=size(temp);
    if si(1)>0
        xnormalvelocity(i+1,1)=temp(1,1);
        xnormalvelocity(i+1,2)=mean(temp(:,3));
        xnormalvelocity(i+1,3)=std(temp(:,3));
        xnormalvelocity(i+1,4)=abs(std(temp(:,3))/mean(temp(:,3)));
    end
end

ynormalvelocity=[];
for i=0:47
    temp=[];
    for j=1:b(1)
        if and(binvelocity(j,2)<(i+1)*10, binvelocity(j,2)>i*10)
            temp=[temp;binvelocity(j,:)];
        end
    end
    si=size(temp);
    if si(1)>0
        ynormalvelocity(i+1,1)=temp(1,2);
        ynormalvelocity(i+1,2)=mean(temp(:,4));
        ynormalvelocity(i+1,3)=std(temp(:,4));
        ynormalvelocity(i+1,4)=abs(std(temp(:,3))/mean(temp(:,3)));
    end
end
% 
% % endpoint=[];
% % endpoint(:,1)=binvelocity(:,1)+binvelocity(:,3);
% % endpoint(:,2)=binvelocity(:,2)+binvelocity(:,4);
% % looplength=size(binvelocity)
% %  figure
% %  for m=1:looplength
% %          vectarrow(binvelocity(m,1:2),endpoint(m,1:2));
% %      hold on;
% %  end
% 
% r=size(xnormalvelocity);
%  figure
%   for t=1:r(1)
%      plot(xnormalvelocity(t,1),xnormalvelocity(t,2),'k')
%      hold on;
%   end
%   
%   s=size(ynormalvelocity);
%    figure
%   for t=1:s(1)
%      plot(ynormalvelocity(t,1),ynormalvelocity(t,2),'k')
%      hold on;
%   end
  
  
  
% 
%  r=size(xnormalvelocity);
%  figure
%   for t=1:r(1)
%      plot(xnormalvelocity(t,1),xnormalvelocity(t,2),'k')
%      hold on;
%   end
%  
%  r=size(ynormalvelocity);
%  figure
%   for t=1:r(1)
%      plot(ynormalvelocity(t,1),ynormalvelocity(t,2),'k')
%      hold on;
%   end
  
  
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



% for t=1:a(1)
%     plot(xdf(t,1),xdf(t,2),'k')
%     plot(ydf(t,1),ydf(t,2),'r')
%     hold on;
% end


%  a=size(velocity);
%  xcenter=velocity(1,1);
%  ycenter=velocity(1,2);
%  xinitialerror=abs(velocity(1,3)-0);
%  yinitialerror=abs(velocity(1,4)-0);
%  for u=2:a(1)
%      xerror=abs(velocity(u,3)-0);
%      yerror=abs(velocity(u,4)-0);
%      if xerror<xinitialerror
%          xcenter=velocity(u,1);
%          xinitialerror=xerror;
%      else 
%          xcenter=xcenter;   
%      end
%      if yerror<yinitialerror
%          ycenter=velocity(u,2);
%          yinitialerror=yerror;
%      else 
%          ycenter=ycenter;   
%      end
%  end
%  
%  velocity(:,1)=abs(velocity(:,1)-xcenter);
%  velocity(:,2)=abs(velocity(:,2)-ycenter);
%  velocity(:,3)=abs(velocity(:,3));
%  velocity(:,4)=abs(velocity(:,4));
%  
%  
%  
%  figure
%   for t=1:a(1)
%      plot(velocity(t,1),velocity(t,3),'k')
%      plot(velocity(t,2),velocity(t,4),'r')
%      hold on;
%  end

% 
% r=size(binvelocity);
% qi=[];
%   for i=1:64
%       qq=[];
%       for t=1:r(1)
%           if and(binvelocity(t,1)<i+1,binvelocity(t,1)>i)
%              qq=[qq;binvelocity(t,:)];
%           end
%       end
%       uu=size(qq);
%       if uu(1)>0
%           aa=max(qq(:,3))-min(qq(:,3));
%           qi(i,1)=(2*i+1)/2;
%           qi(i,2)=aa;
%       end
%   end
% 
%  figure
%  d=size(qi);
%   for t=1:d(1)
%      plot(qi(t,1),qi(t,2),'k')
%      hold on;
%  end

% for i=0:23
%     temp=[];
%     for j=1:b(1)
%         if and(binvelocity(j,2)<i+1, binvelocity(j,2)>i)
%             temp=[temp;binvelocity(j,:)];
%         end
%     end
%     t=size(temp);
% 
%     avgxstd(i+1,1)=temp(1,2);
%     avgxstd(i+1,2)=mean(temp(:,6));
% end
