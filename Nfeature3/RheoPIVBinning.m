num = xlsread('D:\Rheo-PTV\Rheo-PTV Summary\0.3M CTAB 0.3M NaNO3 in D2O\r01.xls','combined');
a=size(num);
data=[];
m=1;
for g=1:a(1)
     xv=[];
     for h=1:a(1)
         if  and(num(h,1)>=(g-1)*0.02,num(h,1)<=g*0.02)
             xv=[xv;num(h,:)];
         end
     end
     qq=size(xv);
     if qq(1)>0
         p=mean(xv(:,2));
         r=std(xv(:,2));
         data(m,1)=(2*g-1)/2*0.02;
         data(m,2)=p;
         data(m,3)=r;
         data(m,4)=qq(1);
         m=m+1;
     end
end