
%pull data from file and read particle positions with Nfeature3
r1=[]
for i=1:254
    im=imread(strcat('C:\Users\Patrick\16-3-10\0.5mlm 60fps',num2str(i-1,'%04.0f'),'.tif'));
    r=Nfeature3(im,1,10,2,15); %(file,lengthscale of noise,size of feature,2,minimum intensity)
%     r(:,2)=480-r(:,2);
    temp=r;
    temp(:,6)=i;
    r1=[r1;temp];
end


s1=r1(:,1:2);
s1(:,3)=r1(:,6);

[lub1]=trackmem(s1,8,2,5,5);%(xyzs,set blocksize,2D,frames feature appears,frames between reappear)

lub=[];
lub=[lub1];

velocity=[];
df=[];
trajectorygap=2;
framerate=60;
for i=1:300
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
binsize=2;
for i=0:512
    for k=0:272
        bin=[];
        for j=1:a(1)
            if and(velocity(j,1)>=i*binsize,velocity(j,1)<=(i+1)*binsize)&& and(velocity(j,2)>=k*binsize,velocity(j,2)<=(k+1)*binsize)
                    bin=[bin;velocity(j,:)];
            end
        end
        y=size(bin);
        if y(1)>0
            binvelocity(t,1)=(2*i+1)/2*binsize;
            binvelocity(t,2)=(2*k+1)/2*binsize;
            binvelocity(t,3)=mean(bin(:,3))*framerate/trajectorygap;
            binvelocity(t,4)=mean(bin(:,4))*framerate/trajectorygap;

%             if binvelocity(t,6)>=100
%                 binvelocity(t,6)=100;
%             end
            t=t+1;
        end
    end
end


mmperpix=0.004;

 s=1;
 x=[];
 y=[];
 c=[];
 x=binvelocity(:,1)*mmperpix;
 y=binvelocity(:,2)*mmperpix;
 c=sqrt((binvelocity(:,3).^2+binvelocity(:,4).^2))*mmperpix;
 
 colormap jet;
 scatter(x,y,s,c,'fill'),colorbar;axis equal;axis([0 1024*mmperpix 0 544*mmperpix]);



        

b=size(binvelocity);

xgradient=[];
m=1;
for g=0:512
    xv=[];
    for h=1:b(1)
        if binvelocity(h,2)==(2*g+1)/2*binsize
            xv=[xv;binvelocity(h,:)];
        end
    end
    uu=size(xv);
    if uu(1)>1
       p=polyfit(xv(:,1),xv(:,3),1);
       r=corrcoef(xv(:,1),xv(:,3));
       xgradient(m,1)=(2*g+1)/2*binsize;
       xgradient(m,2)=p(1);%slope
       xgradient(m,3)=p(2);%intercept
       xgradient(m,4)=abs(r(1,2));
       m=m+1;
    end
end

ygradient=[];
m=1;
for g=0:512
    yv=[];
    for h=1:b(1)
        if binvelocity(h,1)==(2*g+1)/2*binsize
            yv=[yv;binvelocity(h,:)];
        end
    end
    qq=size(yv);
    if qq(1)>1
       p=polyfit(yv(:,2),yv(:,4),1);
       r=corrcoef(yv(:,2),yv(:,4));
       ygradient(m,1)=(2*g+1)/2*binsize;
       ygradient(m,2)=p(1);%slope
       ygradient(m,3)=p(2);%intercept
       ygradient(m,4)=abs(r(1,2));
       
    else
        continue
    end
    m=m+1;
end


