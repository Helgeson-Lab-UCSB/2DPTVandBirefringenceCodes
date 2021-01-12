

b=size(binvelocity);

xgradient=[];
uu=[];

for g=0:95
    xv=[];
    for h=1:b(1)
        if binvelocity(h,2)==(2*g+1)/2*5
            xv=[xv;binvelocity(h,:)];
        end
    end
    p=polyfit(xv(:,1),xv(:,3),1);
    r=corrcoef(xv(:,1),xv(:,3));
    xv(:,3)=p(1);
    xv(:,4)=abs(r(1,2));
    uu=[uu;xv];
end

vv=[];
m=1;
for g=0:127
    yv=[];
    for h=1:b(1)
        if binvelocity(h,1)==(2*g+1)/2*5
            yv=[yv;binvelocity(h,:)];
        end
    end
    qq=size(yv);
    if qq(1)>0
       p=polyfit(yv(:,2),yv(:,4),1);
       r=corrcoef(yv(:,2),yv(:,4));
       yv(:,3)=p(1);
       yv(:,4)=abs(r(1,2));
       vv=[vv;yv];
    end
end


s=10;
x=[];
y=[];
c=[];
x=uu(:,1);
y=uu(:,2);
c=uu(:,3);

colormap jet;
scatter(x,y,s,c,'fill'),colorbar;caxis ([0 0.8]);axis equal;axis([0 640 0 480]);

r=size(uu);
s=size(vv);
gra=[];
m=1;
for t=1:r(1)
    for j=1:s(1)
        if and (uu(t,1)==vv(j,1),uu(t,2)==vv(j,2))
            gra(m,1)=uu(t,1);
            gra(m,2)=uu(t,2);
            gra(m,3)=sqrt((uu(t,3))^2+(vv(j,3))^2);
            m=m+1;
        end
    end
end
            
           
