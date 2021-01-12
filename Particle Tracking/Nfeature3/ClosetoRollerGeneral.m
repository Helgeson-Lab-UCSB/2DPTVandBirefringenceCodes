r1=[]
for i=1:477
    im=imread(strcat('C:\Users\Peng Cheng\Desktop\Measurement1\Measurement10',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,55);
    r(:,2)=480-r(:,2);
    temp=r;
    temp(:,6)=i;
    r1=[r1;temp];
end

% r2=[]
% for i=1:475
%     im=imread(strcat('C:\Users\Peng Cheng\Desktop\Measurement2\Measurement20',num2str(i-1,'%03.0f'),'.tif'));
%     r=Nfeature3(im,1,3,2,40);
%     r(:,2)=480-r(:,2);
%     temp=r;
%     temp(:,6)=i;
%     r2=[r2;temp];
% end
% 
% r3=[]
% for i=1:475
%     im=imread(strcat('C:\Users\Peng Cheng\Desktop\Measurement3\Measurement30',num2str(i-1,'%03.0f'),'.tif'));
%     r=Nfeature3(im,1,3,2,40);
%     r(:,2)=480-r(:,2);
%     temp=r;
%     temp(:,6)=i;
%     r3=[r3;temp];
% end
% 
% r4=[]
% for i=1:475
%     im=imread(strcat('C:\Users\Peng Cheng\Desktop\Measurement4\Measurement40',num2str(i-1,'%03.0f'),'.tif'));
%     r=Nfeature3(im,1,3,2,40);
%     r(:,2)=480-r(:,2);
%     temp=r;
%     temp(:,6)=i;
%     r4=[r4;temp];
% end



s1=r1(:,1:2);
s1(:,3)=r1(:,6);

% s2=r2(:,1:2);
% s2(:,3)=r2(:,6);
% 
% s3=r3(:,1:2);
% s3(:,3)=r3(:,6);
% 
% s4=r4(:,1:2);
% s4(:,3)=r4(:,6);


[lub1]=trackmem(s1,6,2,5,5);
% [lub2]=trackmem(s2,5,2,5,5);
% [lub3]=trackmem(s3,5,2,5,5);
% [lub4]=trackmem(s4,5,2,5,5);
% 
% length1=max(lub1(:,4));
% length2=max(lub2(:,4));
% length3=max(lub3(:,4));
% length4=max(lub4(:,4));
% 
% lub2(:,4)=lub2(:,4)+length1;
% lub3(:,4)=lub3(:,4)+length1+length2;
% lub4(:,4)=lub4(:,4)+length1+length2+length3;
    
lub=[];
lub=[lub1];
% lub=[lub1;lub2;lub3;lub4];


velocity=[];
df=[];
trajectorygap=2;
framerate=60;
for i=1:477
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

velocity(:,1)=velocity(:,1)-612;
velocity(:,2)=velocity(:,2)-466;
len=size(velocity);
result=[];
% outlier=[];
% shearrate=[];
% alpha=1;
% truealpha=1;
m=1;
n=1;
for i=1:len
    x1=velocity(i,1);
    y1=velocity(i,2);
    vx=velocity(i,3)*framerate/trajectorygap;
    vy=velocity(i,4)*framerate/trajectorygap;
    result(m,1)=sqrt(x1^2+y1^2);
    result(m,2)=sqrt(vx^2+vy^2);
%     gamma=(vx-vy)/(x1+y1);
%     alpha=(vx+vy)/(gamma*(x1-y1));
%     if and(abs(gamma)<=0.61,abs(alpha)<=2)  && and(gamma>0,alpha>0)
%        result(m,1)=x1;
%        result(m,2)=y1;
%        result(m,3)=gamma;
%        result(m,4)=alpha;
    m=m+1;
%     else
%        outlier(n,1)=x1;
%        outlier(n,2)=y1;
%        n=n+1;
%     end
%     calc1=60/sqrt(alpha)*asinh(2*sqrt(alpha)*(x2*y1-x1*y2)/((x1^2+y1^2)*(1-alpha)+2*(1+alpha)*x1*y1));
%     calc2=60/sqrt(truealpha)*asinh(2*sqrt(truealpha)*(x2*y1-x1*y2)/((x1^2+y1^2)*(1-truealpha)+2*(1+truealpha)*x1*y1));
%     if and(abs(calc1)<=10,abs(calc2)<=10)
%         shearrate(m,1)=x1;
%         shearrate(m,2)=y1;
%         shearrate(m,3)=calc1;
%         shearrate(m,4)=calc2;
%         m=m+1;
%     end
end

% allresult=[];
% allresult=[result];
% 
% 
% a=size(allresult);
% binvelocity=[];
% t=1;
% for i=0:1500
%     bin=[];
%     for j=1:a(1)
%         if and(allresult(j,1)>=i,allresult(j,1)<=(i+1))
%            bin=[bin;allresult(j,:)];
%         end
%     end
%     y=size(bin);
%     kk=[];
%     if y(1)>0
%        binvelocity(t,1)=(2*i+1)/2;
%        binvelocity(t,2)=mean(bin(:,2));
%        binvelocity(t,3)=std(bin(:,2));
%        kk=size(bin);
%        binvelocity(t,4)=kk(1);
% 
% 
% %             if binvelocity(t,6)>=100
% %                 binvelocity(t,6)=100;
% %             end
%        t=t+1;
%     end
%     
% end

