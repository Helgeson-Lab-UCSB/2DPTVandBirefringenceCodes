r1=[]
for i=1:478
    im=imread(strcat('C:\Users\Peng Cheng\Desktop\Measurement1\Measurement10',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,40);
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
trajectorygap=2;
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


ximagecenter=323;
yimagecenter=236;

velocity(:,1)=velocity(:,1)-ximagecenter;
velocity(:,2)=velocity(:,2)-(480-yimagecenter);


alpha=0.932555;
a1=(-(1+alpha)+2*sqrt(alpha))/(1-alpha);
a2=(-(1+alpha)-2*sqrt(alpha))/(1-alpha);

a1projectx=1/sqrt(1+(a1)^2);
a1projecty=a1/sqrt(1+(a1)^2);
a2projectx=-1/sqrt(1+(a2)^2);
a2projecty=-a2/sqrt(1+(a2)^2);

vv=[];
vv(:,1)=(velocity(:,2)-a2*velocity(:,1))/(a1-a2)*sqrt(a1^2+1);
vv(:,2)=-(velocity(:,2)-a1*velocity(:,1))/(a2-a1)*sqrt(a2^2+1);
vv(:,3)=velocity(:,3)*a1projectx+velocity(:,4)*a1projecty;
vv(:,4)=velocity(:,3)*a2projectx+velocity(:,4)*a2projecty;



a=size(vv);

xlowlim=round(min(vv(:,1)));
xhighlim=round(max(vv(:,1)));
ylowlim=round(min(vv(:,2)));
yhighlim=round(max(vv(:,2)));

xbin=round((xhighlim-xlowlim)/5);
ybin=round((yhighlim-ylowlim)/5);

binvelocity=[];
t=1;
for i=0:xbin
    for k=0:ybin
        bin=[];
        for j=1:a(1)
            if and(vv(j,1)>=i*5+xlowlim,vv(j,1)<=(i+1)*5+xlowlim)&& and(vv(j,2)>=k*5+ylowlim,vv(j,2)<=(k+1)*5+ylowlim)
                    bin=[bin;vv(j,:)];
            end
        end
        y=size(bin);
        if y(1)>0
            binvelocity(t,1)=(2*i+1)/2*5+xlowlim;
            binvelocity(t,2)=(2*k+1)/2*5+ylowlim;
            binvelocity(t,3)=mean(bin(:,3))*framerate/trajectorygap;
            binvelocity(t,4)=mean(bin(:,4))*framerate/trajectorygap;

%             if binvelocity(t,6)>=100
%                 binvelocity(t,6)=100;
%             end
            t=t+1;
        end
    end
end

b=size(binvelocity);

xgradient=[];
m=1;
for g=0:300
    xv=[];
    for h=1:b(1)
        if binvelocity(h,2)==(2*g+1)/2*5+ylowlim
            xv=[xv;binvelocity(h,:)];
        else
            continue
        end
    end
    qq=size(xv);
    if qq(1)>1
       p=polyfit(xv(:,1),xv(:,3),1);
       r=corrcoef(xv(:,1),xv(:,3));
       xgradient(m,1)=(2*g+1)/2*5+ylowlim;
       xgradient(m,2)=p(1);
       xgradient(m,3)=p(2);
       xgradient(m,4)=abs(r(1,2));
       m=m+1;
    else 
        continue
    end
end

ygradient=[];
m=1;
for g=0:300
    yv=[];
    for h=1:b(1)
        if binvelocity(h,1)==(2*g+1)/2*5+xlowlim;
            yv=[yv;binvelocity(h,:)];
        else
            continue
        end
    end
    qq=size(yv);
    if qq(1)>0
       p=polyfit(yv(:,2),yv(:,4),1);
       r=corrcoef(yv(:,2),yv(:,4));
       ygradient(m,1)=(2*g+1)/2*5+xlowlim;
       ygradient(m,2)=p(1);
       ygradient(m,3)=p(2);
       ygradient(m,4)=abs(r(1,2));
       m=m+1;
    else 
        continue
    end
end



