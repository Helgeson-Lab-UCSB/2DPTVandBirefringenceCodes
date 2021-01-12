r1=[]
for i=1:475
    im=imread(strcat('C:\Users\Peng Cheng\Desktop\Measurement1\Measurement10',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,60);
    r(:,2)=480-r(:,2);
    temp=r;
    temp(:,6)=i;
    r1=[r1;temp];
end

% r2=[]
% for i=1:445
%     im=imread(strcat('C:\Users\HelgesonGp\Desktop\Measurement2\Measurement20',num2str(i-1,'%03.0f'),'.tif'));
%     r=Nfeature3(im,1,3,2,65);
%     r(:,2)=640-r(:,2);
%     temp=r;
%     temp(:,6)=i;
%     r2=[r2;temp];
% end

% r3=[]
% for i=1:117
%     im=imread(strcat('C:\Users\HelgesonGp\Desktop\Measurement3\Measurement30',num2str(i-1,'%03.0f'),'.tif'));
%     r=Nfeature3(im,1,3,2,65);
%     r(:,2)=640-r(:,2);
%     temp=r;
%     temp(:,6)=i;
%     r3=[r3;temp];
% end
% 
% r4=[]
% for i=1:117
%     im=imread(strcat('C:\Users\HelgesonGp\Desktop\Measurement4\Measurement40',num2str(i-1,'%03.0f'),'.tif'));
%     r=Nfeature3(im,1,3,2,65);
%     r(:,2)=640-r(:,2);
%     temp=r;
%     temp(:,6)=i;
%     r4=[r4;temp];
% end



s1=r1(:,1:2);
s1(:,3)=r1(:,6);

% s2=r2(:,1:2);
% s2(:,3)=r2(:,6);

% s3=r3(:,1:2);
% s3(:,3)=r3(:,6);
% 
% s4=r4(:,1:2);
% s4(:,3)=r4(:,6);


[lub1]=trackmem(s1,5,2,5,2);
% [lub2]=trackmem(s2,5,2,5,5);
% [lub3]=trackmem(s3,5,2,5,5);
% [lub4]=trackmem(s4,5,2,5,5);

% length1=max(lub1(:,4));
% length2=max(lub2(:,4));
% length3=max(lub3(:,4));
% length4=max(lub4(:,4));

% lub2(:,4)=lub2(:,4)+length1;
% lub3(:,4)=lub3(:,4)+length1+length2;
% lub4(:,4)=lub4(:,4)+length1+length2+length3;
    
lub=[];
lub=[lub1];
% lub=[lub1;lub2;lub3;lub4];

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

particlenumber=max(lub(:,4));
good=[];
m=1;
framerate=60;
for i=1:particlenumber
    goodparticle=lub(find(lub(:,4)==i),:);
    if range(goodparticle(:,3))>=60
       start1(:,4)=m;
       good=[good;goodparticle];
       m=m+1;
    end          
end

ximagecenter=313;
yimagecenter=237;
gcalc=0.123;
good(:,1)=good(:,1)-ximagecenter;
good(:,2)=good(:,2)-(480-yimagecenter);

goodparticle=max(good(:,4));
trajectorygap=20;
correctparticle=[];
for i=1:goodparticle
    p1=good(find(good(:,4)==i),:);
    a=size(p1);
    temp=[];
    m=1;
    alpha=1;
    for k=1:a(1)-trajectorygap
        start=p1(k,:);
        finish=p1(k+trajectorygap,:);
        x1=start(1,1);
        y1=start(1,2);
        x2=(finish(1,1)-start(1,1))/(finish(1,3)-start(1,3));
        y2=(finish(1,2)-start(1,2))/(finish(1,3)-start(1,3));
        calc1=framerate/(finish(1,3)-start(1,3))/sqrt(alpha)*asinh(2*sqrt(alpha)*(x2*y1-x1*y2)/((x1^2+y1^2)*(1-alpha)+2*(1+alpha)*x1*y1));
        temp(m,1)=calc1;
        m=m+1;
    end
    kkkk=size(temp);
    if (kkkk(1)>0) && and(mean(temp(:,1))>=0,mean(temp(:,1))<=2*gcalc)
        correctparticle=[correctparticle;p1];
    end
end

u=max(correctparticle(:,4));
firstpeak=[];
secondpeak=[];
dataplot=[];
m=1;
for i=1:u
    p1=correctparticle(find(correctparticle(:,4)==i),:);
    r=size(p1);
    g=[];
    m=1;
    alpha=1;
    for j=1:r(1)-1
      x1=p1(j,1);
      y1=p1(j,2);
      vx=(p1(j+1,1)-p1(j,1))/(p1(j+1,3)-p1(j,3));
      vy=(p1(j+1,2)-p1(j,2))/(p1(j+1,3)-p1(j,3));
      calc1=(x1*vx-y1*vy)/sqrt(x1^2+y1^2);
      g(m,1)=p1(j,3);
      g(m,2)=calc1;
      m=m+1;
    end
    llll=size(g);
    if llll(1)>0
       x=fft(g(:,2));
       x(1)=[];  
       n=length(x);
       power = abs(x(1:floor(n/2))).^2;
       nyquist=framerate/2;
       freq=(1:n/2)/(n/2)*nyquist;
       spectrum=[];
       spectrum(:,1)=freq';
       spectrum(:,2)=power;
       t=size(spectrum);
       maximum=0;
       temp=[];
       for k=1:t(1)
           if spectrum(k,2)>maximum
              maximum=spectrum(k,2);
              temp=spectrum(k,:);
              maxloc=spectrum(k,1);
           end
       end
       firstpeak=[firstpeak;temp];
       p1(:,3)=maximum;
       p1(:,4)=maxloc;
       dataplot=[dataplot;p1];
%        secondmaximum=0;
%        temp2=[];
%        for o=1:t(1)
%            if (abs(o-maxloc)>=3) && and(spectrum(o,2)>secondmaximum,spectrum(o,2)<maximum)
%                secondmaximum=spectrum(o,2);
%                temp2=spectrum(o,:);
%            end
%        end
%        secondpeak=[secondpeak;temp2];
    end

end  
    
    
hu=size(firstpeak);
avg=mean(firstpeak(:,1));
cleanfirstpeak=[];
m=1;
for i=1:hu(1)
    if and(firstpeak(i,1)<=1.5*avg,firstpeak(i,1)>=0.5*avg)
        cleanfirstpeak(m,:)=firstpeak(i,:);
        m=m+1;
    end
end


% ou=size(secondpeak);
% avg2=mean(secondpeak(:,1));
% cleansecondpeak=[];
% m=1;
% for i=1:ou(1)
%     if secondpeak(i,1)<=1.5*avg2
%         cleansecondpeak(m,:)=secondpeak(i,:);
%         m=m+1;
%     end  
% end


mean(cleanfirstpeak(:,1))

std(cleanfirstpeak(:,1))

mean(cleanfirstpeak(:,2))

std(cleanfirstpeak(:,2))

a=size(firstpeak);
b=size(cleanfirstpeak);
b(1)/a(1)
    
    
% mean(cleansecondpeak(:,1))
% 
% std(cleansecondpeak(:,1))
% 
% mean(cleansecondpeak(:,2))
% 
% std(cleansecondpeak(:,2))    

% s=1;
% x=[];
% y=[];
% c=[];
% x=dataplot(:,1)+ximagecenter;
% y=dataplot(:,2)+(480-yimagecenter);
% c=dataplot(:,3);
% 
% colormap jet;
% scatter(x,y,s,c,'fill'),colorbar;axis equal;axis([0 640 0 480]);

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
 
%  velocity(:,1)=abs(velocity(:,1)-xcenter);
%  velocity(:,2)=abs(velocity(:,2)-ycenter);
%  velocity(:,3)=abs(velocity(:,3));
%  velocity(:,4)=abs(velocity(:,4));


% a=size(velocity);
% 
% % binvelocity=[];
% % t=1;
% % for i=0:63
% %     for k=0:47
% %         bin=[];
% %         for j=1:a(1)
% %             if and(velocity(j,1)>=i*10,velocity(j,1)<=(i+1)*10)&& and(velocity(j,2)>=k*10,velocity(j,2)<=(k+1)*10)
% %                     bin=[bin;velocity(j,:)];
% %             end
% %         end
% %         y=size(bin);
% %         if y(1)>0
% %             binvelocity(t,1)=(2*i+1)/2;
% %             binvelocity(t,2)=(2*k+1)/2;
% %             binvelocity(t,3)=mean(bin(:,3));
% %             binvelocity(t,4)=mean(bin(:,4));
% % 
% % %             if binvelocity(t,6)>=100
% % %                 binvelocity(t,6)=100;
% % %             end
% %             t=t+1;
% %         end
% %     end
% % end
% 
% % b=size(binvelocity);
% % xnormalvelocity=[];
% % for i=0:63
% %     temp=[];
% %     for j=1:b(1)
% %         if and(binvelocity(j,1)<i+1, binvelocity(j,1)>i)
% %             temp=[temp;binvelocity(j,:)];
% %         end
% %     end
% %     si=size(temp);
% %     if si(1)>0
% %         xnormalvelocity(i+1,1)=temp(1,1);
% %         xnormalvelocity(i+1,2)=mean(temp(:,3));
% %         xnormalvelocity(i+1,3)=std(temp(:,3));
% %         xnormalvelocity(i+1,4)=abs(std(temp(:,3))/mean(temp(:,3)));
% %     end
% % end
% % 
% % ynormalvelocity=[];
% % for i=0:47
% %     temp=[];
% %     for j=1:b(1)
% %         if and(binvelocity(j,2)<i+1, binvelocity(j,2)>i)
% %             temp=[temp;binvelocity(j,:)];
% %         end
% %     end
% %     si=size(temp);
% %     if si(1)>0
% %         ynormalvelocity(i+1,1)=temp(1,2);
% %         ynormalvelocity(i+1,2)=mean(temp(:,4));
% %         ynormalvelocity(i+1,3)=std(temp(:,4));
% %         ynormalvelocity(i+1,4)=abs(std(temp(:,3))/mean(temp(:,3)));
% %     end
% % end
% % 
% % % endpoint=[];
% % % endpoint(:,1)=binvelocity(:,1)+binvelocity(:,3);
% % % endpoint(:,2)=binvelocity(:,2)+binvelocity(:,4);
% % % looplength=size(binvelocity)
% % %  figure
% % %  for m=1:looplength
% % %          vectarrow(binvelocity(m,1:2),endpoint(m,1:2));
% % %      hold on;
% % %  end
% % 
% 
% ximagecenter=319;
% yimagecenter=248;
% velocity(:,1)=velocity(:,1)-ximagecenter;
% velocity(:,2)=velocity(:,2)-(480-yimagecenter);
% len=max(velocity(:,5));
% result=[];
% % shearrate=[];
% % alpha=1;
% % truealpha=1;
% m=1;
% n=1;
% for i=1:len
%     
%     transientvelocity=velocity(find(velocity(:,5)==i),:);
%     a=size(transientvelocity);
%     transientresult=[];
%     m=1;
%     for j=1:a
%         x1=transientvelocity(j,1);
%         y1=transientvelocity(j,2);
%         vx=transientvelocity(j,3)*framerate/trajectorygap;
%         vy=transientvelocity(j,4)*framerate/trajectorygap;
%         gamma=(vx-vy)/(x1+y1);
%         alpha=(vx+vy)/(gamma*(x1-y1));
%         if and(abs(gamma)<=0.08,abs(alpha)<=2)&& alpha>=0
%            transientresult(m,1)=x1;
%            transientresult(m,2)=y1;
%            transientresult(m,3)=gamma;
%            transientresult(m,4)=alpha;
%            
%            m=m+1;
%         end
%     end
%     
%     result(i,1)=mean(transientresult(:,3));
%     result(i,2)=std(transientresult(:,3));
%     result(i,3)=mean(transientresult(:,4));
%     result(i,4)=std(transientresult(:,4));
%     result(i,5)=i/framerate;
%     result(i,6)=m/a(1); %acceptance ratio
% %     calc1=60/sqrt(alpha)*asinh(2*sqrt(alpha)*(x2*y1-x1*y2)/((x1^2+y1^2)*(1-alpha)+2*(1+alpha)*x1*y1));
% %     calc2=60/sqrt(truealpha)*asinh(2*sqrt(truealpha)*(x2*y1-x1*y2)/((x1^2+y1^2)*(1-truealpha)+2*(1+truealpha)*x1*y1));
% %     if and(abs(calc1)<=10,abs(calc2)<=10)
% %         shearrate(m,1)=x1;
% %         shearrate(m,2)=y1;
% %         shearrate(m,3)=calc1;
% %         shearrate(m,4)=calc2;
% %         m=m+1;
% %     end
% end


% r=size(outlier);
% im=imread('C:\Users\HelgesonGp\Desktop\03-15-2013 Strong Flow Calibration\Eta=0.4.tif');  
% figure
% imagesc(im);colormap gray;axis equal
% hold on;
%    for t=1:r(1)
%       plot((outlier(t,1)+ximagecenter),(640-(outlier(t,2)+(640-yimagecenter))),'r')
%       
%       hold on;
%    end
% r=size(binvelocity);
%  figure
%   for t=1:r(1)
%      plot(binvelocity(t,1)*10,binvelocity(t,3)*60,'k')
%      plot(binvelocity(t,2)*10,binvelocity(t,4)*60,'r')
%      hold on;
%   end
%   
  

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
% m=1;
% bigerror=[];
% len=size(result);
% for i=1:len(1)
%     error=abs(result(i,3)-mean(result(:,3)));
%     if error>=std(result(:,3))
%         bigerror(m,1)=result(i,1);
%         bigerror(m,2)=result(i,2);
%         m=m+1;
%     end
% end