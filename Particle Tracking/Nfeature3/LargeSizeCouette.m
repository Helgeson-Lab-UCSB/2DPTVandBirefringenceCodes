
r1=[]
for i=1:99
    im=imread(strcat('C:\Users\Peng Cheng\Desktop\Measurement1\Measurement10',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,15);
    r(:,2)=640-r(:,2);
    temp=r;
    temp(:,6)=i;
    r1=[r1;temp];
end 


s1=r1(:,1:2);
s1(:,3)=r1(:,6);



[lub1]=trackmem(s1,5,2,5,5);
lub=[];
lub=[lub1];

velocity=[];
df=[];
trajectorygap=1;
framerate=20;
for i=1:99
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

velocity(:,1)=velocity(:,1)+1.65*286-308;
velocity(:,2)=velocity(:,2)+1.65*286-373;
len=size(velocity);
result3=[];

m=1;
n=1;
for i=1:len
    x1=velocity(i,1);
    y1=velocity(i,2);
    vx=velocity(i,3)*framerate/trajectorygap;
    vy=velocity(i,4)*framerate/trajectorygap;
    result3(m,1)=sqrt(x1^2+y1^2);
    result3(m,2)=sqrt(vx^2+vy^2);

    m=m+1;

end

allresult=[];
allresult=[result3];


a=size(allresult);
binvelocity=[];
t=1;
for i=0:600
    bin=[];
    for j=1:a(1)
        if and(allresult(j,1)>=i,allresult(j,1)<=(i+1))
           bin=[bin;allresult(j,:)];
        end
    end
    y=size(bin);
    kk=[];
    if y(1)>0
       binvelocity(t,1)=(2*i+1)/2;
       binvelocity(t,2)=mean(bin(:,2));
       binvelocity(t,3)=std(bin(:,2));
       kk=size(bin);
       binvelocity(t,4)=kk(1);
       t=t+1;
    end
    
end