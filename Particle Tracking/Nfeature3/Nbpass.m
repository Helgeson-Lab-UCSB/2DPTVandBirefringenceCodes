function s = Nbpass(img,lambda,w)
clear s t
a=double(img);
[xdim ydim]=size(a);
N=xdim*ydim;
anew=zeros(N,1);
for y=1:ydim
    for x=1:xdim
        m=xdim*(y-1)+x;
        anew(m)=a(x,y);
    end
end
% xdim=length(a(:,1));
% ydim=length(a(1,:));
pref=1/(2*w+1)^2.0;
B=sum(exp(-([-w:w].^2/(4*lambda^2))))^2;  %normalization
Aw=zeros(N,1);
Al=zeros(N,1);
% Aw=zeros(xdim,ydim);
% Al=zeros(xdim,ydim);
for x=1+w:xdim-w,
    for y=1+w:ydim-w,
        for i=-w:w,
            for j=-w:w,
        Aw(xdim*(y-1)+x,1)=Aw(xdim*(y-1)+x,1)+anew(xdim*(y+j-1)+x+i,1);
        Al(xdim*(y-1)+x,1)=Al(xdim*(y-1)+x,1)+anew(xdim*(y+j-1)+x+i,1)*exp(-(i^2+j^2)/(4*lambda^2));
%         for i=-w:w,
%             for j=-w:w,
%                 Aw(x,y)=Aw(x,y)+a(x+i,y+j);  %backgound without prefactor
%                 Al(x,y)=Al(x,y)+a(x+i,y+j)*exp(-(i^2+j^2)/(4*lambda^2));  %convolusion with gaussian surface
%             %without devided by B    
            end
        end
    end
end
bg=Aw*pref;  %background
nr=Al/B;     
s=nr-bg;   %difference between noise reduced and background images: estimate of real image.
