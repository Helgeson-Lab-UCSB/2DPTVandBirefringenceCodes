function s = bpass(img,lambda,w)
clear s t
a=double(img);
xdim=length(a(:,1));
ydim=length(a(1,:));
pref=1/(2*w+1)^2;
B=sum(exp(-([-w:w].^2/(4*lambda^2))))^2;  %normalization
Aw=zeros(xdim,ydim);
Al=zeros(xdim,ydim);
for x=1+w:xdim-w,
    for y=1+w:ydim-w,
%         Aw(x,y)=sum(sum(a(x-w:x+w,y-w:y+w),1));
%         Al(x,y)=sum(sum(a(x-w:x+w,y-w:y+w).*exp(-((ones(2*w+1,1)*[-w:w]).^2+([-w:w]'*ones(1,2*w+1)).^2)/(4*lambda^2)),1));
        for i=-w:w,
            for j=-w:w,
                Aw(x,y)=Aw(x,y)+a(x+i,y+j);  %backgound without prefactor
                Al(x,y)=Al(x,y)+a(x+i,y+j)*exp(-(i^2+j^2)/(4*lambda^2));  %convolusion with gaussian surface
            %without devided by B    
            end
        end
    end
end
bg=Aw*pref;  %background
nr=Al/B;     
s=nr-bg;   %difference between noise reduced and background images: estimate of real image.
