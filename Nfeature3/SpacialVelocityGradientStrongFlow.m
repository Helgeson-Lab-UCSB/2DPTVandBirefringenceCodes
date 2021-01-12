r1=[]
for i=1:478
    im=imread(strcat('C:\Users\Peng Cheng\Desktop\Measurement1\Measurement10',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,60);
    r(:,2)=480-r(:,2);
    temp=r;
    temp(:,6)=i;
    r1=[r1;temp];
end

s1=r1(:,1:2);
s1(:,3)=r1(:,6);

[lub1]=trackmem(s1,6,2,5,5);

lub=[];
lub=[lub1];

velocity=[];
df=[];
trajectorygap=3;
framerate=30;
for i=1:478
    start1=lub(find(lub(:,3)==i),:);
    finish1=lub(find(lub(:,3)==i+trajectorygap),:);
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

binvelocity=[];
t=1;
for i=0:127
    for k=0:95
        bin=[];
        for j=1:a(1)
            if and(velocity(j,1)>=i*5,velocity(j,1)<=(i+1)*5)&& and(velocity(j,2)>=k*5,velocity(j,2)<=(k+1)*5)
                    bin=[bin;velocity(j,:)];
            end
        end
        y=size(bin);
        if y(1)>0
            binvelocity(t,1)=(2*i+1)/2*5;
            binvelocity(t,2)=(2*k+1)/2*5;
            binvelocity(t,3)=mean(bin(:,3))*framerate/trajectorygap;
            binvelocity(t,4)=mean(bin(:,4))*framerate/trajectorygap;

%             if binvelocity(t,6)>=100
%                 binvelocity(t,6)=100;
%             end
            t=t+1;
        end
    end
end


ximagecenter=317;
yimagecenter=269;

binvelocity(:,1)=binvelocity(:,1)-ximagecenter;
binvelocity(:,2)=binvelocity(:,2)-(480-yimagecenter);



alpha=0.277982;
a1=(-(1+alpha)+2*sqrt(alpha))/(1-alpha);
a2=(-(1+alpha)-2*sqrt(alpha))/(1-alpha);

a1projectx=1/sqrt(1+(a1)^2);
a1projecty=a1/sqrt(1+(a1)^2);
a2projectx=-1/sqrt(1+(a2)^2);
a2projecty=-a2/sqrt(1+(a2)^2);

vv=[];
vv(:,1)=(binvelocity(:,2)-a2*binvelocity(:,1))/(a1-a2)*sqrt(a1^2+1);
vv(:,2)=-(binvelocity(:,2)-a1*binvelocity(:,1))/(a2-a1)*sqrt(a2^2+1);
vv(:,3)=binvelocity(:,3)*a1projectx+binvelocity(:,4)*a1projecty;
vv(:,4)=binvelocity(:,3)*a2projectx+binvelocity(:,4)*a2projecty;

 b=size(vv);
 
xlowlim=round(min(vv(:,1)));
xhighlim=round(max(vv(:,1)));
ylowlim=round(min(vv(:,2)));
yhighlim=round(max(vv(:,2)));


xgradient=[];
m=1;
for g=(ylowlim+20):yhighlim
    
    xv=[];
    for h=1:b(1)
        if  and(vv(h,2)>=g,vv(h,2)<=g+1)
            xv=[xv;vv(h,:)];
        end
    end
    qq=size(xv);
     if qq(1)>0
       p=polyfit(xv(:,1),xv(:,3),1);
       r=corrcoef(xv(:,1),xv(:,3));
       xgradient(m,1)=(2*g+1)/2;
       xgradient(m,2)=p(1);
       xgradient(m,3)=p(2);
       xgradient(m,4)=abs(r(1,2));
       m=m+1;
     end
end

 b=size(vv);

ygradient=[];
m=1;
for g=(xlowlim+30):xhighlim
    
    xv=[];
    for h=1:b(1)
        if  and(vv(h,1)>=g,vv(h,1)<=g+1)
            xv=[xv;vv(h,:)];
        end
    end
    qq=size(xv);
     if qq(1)>0
       p=polyfit(xv(:,2),xv(:,4),1);
       r=corrcoef(xv(:,2),xv(:,4));
       ygradient(m,1)=(2*g+1)/2;
       ygradient(m,2)=p(1);
       ygradient(m,3)=p(2);
       ygradient(m,4)=abs(r(1,2));
       m=m+1;
     end
end


