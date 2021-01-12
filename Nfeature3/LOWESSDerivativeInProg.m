a=size(velocity,1);
fitmat=zeros(a,5);
p=200;
tic
for i=1:1
    g=1;
    for m=1:a
        r=sqrt((velocity(i,1)-velocity(m,1))^2+(velocity(i,2)-velocity(m,2))^2);
        if(r<p)        
            fitmat(g,1)=velocity(m,1);
            fitmat(g,2)=velocity(m,2);
            fitmat(g,3)=velocity(m,3);
            fitmat(g,4)=velocity(m,4);
            fitmat(g,5)=(1-(r/p)^3)^3;
            fitmat(g,6)=r;
            g=g+1;
        end
    end
    fitmat(any(fitmat==0,2),:)=[];
    fu=fit([fitmat(:,1) fitmat(:,2)],fitmat(:,3),'poly11','Weights',fitmat(:,5));
    fv=fit([fitmat(:,1) fitmat(:,2)],fitmat(:,4),'poly11','Weights',fitmat(:,5));
    ucoeff=coeffvalues(fu);
    vcoeff=coeffvalues(fv);
    dpux(i,1)=velocity(i,1);
    dpuy(i,1)=velocity(i,1);
    dpvx(i,1)=velocity(i,1);
    dpvy(i,1)=velocity(i,1);
    dpux(i,2)=velocity(i,2);
    dpuy(i,2)=velocity(i,2);
    dpvx(i,2)=velocity(i,2);
    dpvy(i,2)=velocity(i,2);
    dpux(i,3)=ucoeff(2)*framerate;
    dpuy(i,3)=ucoeff(3)*framerate;
    dpvx(i,3)=vcoeff(2)*framerate;
    dpvy(i,3)=vcoeff(3)*framerate;
end
toc