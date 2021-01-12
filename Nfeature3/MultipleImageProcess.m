r1=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement1\Measurement10',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r1=[r1;temp];
end

r2=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement2\Measurement20',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp(:,6)=i;
    r2=[r2;temp];
end

r3=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement3\Measurement30',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r3=[r3;temp];
end

r4=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement4\Measurement40',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r4=[r4;temp];
end

r5=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement5\Measurement50',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r5=[r5;temp];
end
r6=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement6\Measurement60',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r6=[r6;temp];
end
r7=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement7\Measurement70',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r7=[r7;temp];
end
r8=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement8\Measurement80',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r8=[r8;temp];
end
r9=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement9\Measurement90',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r9=[r9;temp];
end
r10=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement10\Measurement100',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r10=[r10;temp];
end

r11=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement11\Measurement110',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r11=[r11;temp];
end

r12=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement12\Measurement120',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r12=[r12;temp];
end

r13=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement13\Measurement130',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r13=[r13;temp];
end

r14=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement14\Measurement140',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r14=[r14;temp];
end

r15=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement15\Measurement150',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r15=[r15;temp];
end


r16=[]
for i=1:118
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Low Density Particle\Measurement14\Measurement140',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,100);
    temp=r;
    temp(:,6)=i;
    r16=[r16;temp];
end

s1=r1(:,1:2);
s1(:,3)=r1(:,6);

s2=r2(:,1:2);
s2(:,3)=r2(:,6);

s3=r3(:,1:2);
s3(:,3)=r3(:,6);

s4=r4(:,1:2);
s4(:,3)=r4(:,6);

s5=r5(:,1:2);
s5(:,3)=r5(:,6);
s6=r6(:,1:2);
s6(:,3)=r6(:,6);
s7=r7(:,1:2);
s7(:,3)=r7(:,6);
s8=r8(:,1:2);
s8(:,3)=r8(:,6);
s9=r9(:,1:2);
s9(:,3)=r9(:,6);
s10=r10(:,1:2);
s10(:,3)=r10(:,6);


s11=r11(:,1:2);
s11(:,3)=r11(:,6);

s12=r12(:,1:2);
s12(:,3)=r12(:,6);

s13=r13(:,1:2);
s13(:,3)=r13(:,6);

s14=r14(:,1:2);
s14(:,3)=r14(:,6);

s15=r15(:,1:2);
s15(:,3)=r15(:,6);

s16=r16(:,1:2);
s16(:,3)=r16(:,6);

[lub1]=trackmem(s1,5,2,5,5);
[lub2]=trackmem(s2,5,2,5,5);
[lub3]=trackmem(s3,5,2,5,5);
[lub4]=trackmem(s4,5,2,5,5);
[lub5]=trackmem(s5,5,2,5,5);
[lub6]=trackmem(s6,5,2,5,5);
[lub7]=trackmem(s7,5,2,5,5);
[lub8]=trackmem(s8,5,2,5,5);
[lub9]=trackmem(s9,5,2,5,5);
[lub10]=trackmem(s10,5,2,5,5);
[lub11]=trackmem(s11,5,2,5,5);
[lub12]=trackmem(s12,5,2,5,5);
[lub13]=trackmem(s13,5,2,5,5);
[lub14]=trackmem(s14,5,2,5,5);
[lub15]=trackmem(s15,5,2,5,5);
[lub16]=trackmem(s16,5,2,5,5);


length1=max(lub1(:,4));
length2=max(lub2(:,4));
length3=max(lub3(:,4));
length4=max(lub4(:,4));
length5=max(lub5(:,4));
length6=max(lub6(:,4));
length7=max(lub7(:,4));
length8=max(lub8(:,4));
length9=max(lub9(:,4));
length10=max(lub10(:,4));
length11=max(lub11(:,4));
length12=max(lub12(:,4));
length13=max(lub13(:,4));
length14=max(lub14(:,4));
length15=max(lub15(:,4));

lub2(:,4)=lub2(:,4)+length1;
lub3(:,4)=lub3(:,4)+length1+length2;
lub4(:,4)=lub4(:,4)+length1+length2+length3;
lub5(:,4)=lub5(:,4)+length1+length2+length3+length4;
lub6(:,4)=lub6(:,4)+length1+length2+length3+length4+length5;
lub7(:,4)=lub7(:,4)+length1+length2+length3+length4+length5+length6;
lub8(:,4)=lub8(:,4)+length1+length2+length3+length4+length5+length6+length7;
lub9(:,4)=lub9(:,4)+length1+length2+length3+length4+length5+length6+length7+length8;
lub10(:,4)=lub10(:,4)+length1+length2+length3+length4+length5+length6+length7+length8+length9;
lub11(:,4)=lub11(:,4)+length1+length2+length3+length4+length5+length6+length7+length8+length9+length10;
lub12(:,4)=lub12(:,4)+length1+length2+length3+length4+length5+length6+length7+length8+length9+length10+length11;
lub13(:,4)=lub13(:,4)+length1+length2+length3+length4+length5+length6+length7+length8+length9+length10+length11+length12;
lub14(:,4)=lub14(:,4)+length1+length2+length3+length4+length5+length6+length7+length8+length9+length10+length11+length12+length13;
lub15(:,4)=lub15(:,4)+length1+length2+length3+length4+length5+length6+length7+length8+length9+length10+length11+length12+length13+length14;
lub16(:,4)=lub16(:,4)+length1+length2+length3+length4+length5+length6+length7+length8+length9+length10+length11+length12+length13+length14+length15;

lub=[];
lub=[lub1;lub2;lub3;lub4;lub5;lub6;lub7;lub8;lub9;lub10;lub11;lub12;lub13;lub14;lub15;lub16];

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

totallength=size(lub);
for k=1:totallength
    pt=lub(find(lub(:,4)==k),:);
    pt=sortrows(pt,3);
    hold on;
    plot(pt(:,1),pt(:,2),'b')
end
box on;

velocity=[];
df=[];
for i=1:117
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

b=size(binvelocity);
xnormalvelocity=[];
for i=0:63
    temp=[];
    for j=1:b(1)
        if and(binvelocity(j,1)<i+1, binvelocity(j,1)>i)
            temp=[temp;binvelocity(j,:)];
        end
    end
    xnormalvelocity(i+1,1)=temp(1,1);
    xnormalvelocity(i+1,2)=mean(temp(:,3));
    xnormalvelocity(i+1,3)=std(temp(:,3));
    xnormalvelocity(i+1,4)=abs(std(temp(:,3))/mean(temp(:,3)));
end

ynormalvelocity=[];
for i=0:47
    temp=[];
    for j=1:b(1)
        if and(binvelocity(j,2)<i+1, binvelocity(j,2)>i)
            temp=[temp;binvelocity(j,:)];
        end
    end
    ynormalvelocity(i+1,1)=temp(1,2);
    ynormalvelocity(i+1,2)=mean(temp(:,4));
    ynormalvelocity(i+1,3)=std(temp(:,4));
    ynormalvelocity(i+1,4)=abs(std(temp(:,3))/mean(temp(:,3)));
end



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
% % for t=1:a(1)
% %     plot(xdf(t,1),xdf(t,2),'k')
% %     plot(ydf(t,1),ydf(t,2),'r')
% %     hold on;
% % end
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
%  figure
%   for t=1:a(1)
%      plot(velocity(t,1),velocity(t,3),'k')
%      plot(velocity(t,2),velocity(t,4),'r')
%      hold on;
%  end