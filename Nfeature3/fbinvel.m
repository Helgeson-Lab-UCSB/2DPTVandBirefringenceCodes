function [ outputvel, outputvelmatu,outputvelmatv ] = fbinvel(inputvel,imagexres,imageyres,velbinsize,edgecutoff)
%fbinvel bins a velocity matrix
%   Takes a matrix inputvel organized as xpos, ypos, xvel, yvel and returns
%   outputvel (a binned matrix of similar form) and
%   outputvelmat[xpos,ypos]=outputvel(xpos,ypos)
xres=1;
yres=1;

a=size(inputvel);
outputvelmatu=zeros(imagexres/(yres*velbinsize),imageyres/(yres*velbinsize));
outputvelmatv=zeros(imagexres/(yres*velbinsize),imageyres/(yres*velbinsize));
outputvel=zeros(imagexres*imageyres/(velbinsize*velbinsize),3);

tic
counter=zeros(ceil(imagexres/(xres*velbinsize)),ceil(imageyres/(yres*velbinsize)));
t=1;
for j=1:a(1)%pivot data to (x/binsize+1,y/binsize+1) matrix rounding to nearest integer and adding 1
    xclosest=round(inputvel(j,1)/velbinsize)+1;
    yclosest=round(inputvel(j,2)/velbinsize)+1;
    outputvelmatu(xclosest,yclosest)=outputvelmatu(xclosest,yclosest)+inputvel(j,3);
    outputvelmatv(xclosest,yclosest)=outputvelmatv(xclosest,yclosest)+inputvel(j,4);
    counter(xclosest,yclosest)=counter(xclosest,yclosest)+1;
end
for i=edgecutoff:ceil(imagexres/velbinsize)-edgecutoff %average values in matrix elements
    for k=edgecutoff:ceil(imageyres/velbinsize)-edgecutoff
        if counter(i+1,k+1)>0
            outputvel(t,1)=i*velbinsize;
            outputvel(t,2)=k*velbinsize;
            outputvel(t,3)=outputvelmatu(i+1,k+1)/counter(i+1,k+1);
            outputvel(t,4)=outputvelmatv(i+1,k+1)/counter(i+1,k+1);
            outputvelmatu(i+1,k+1)=outputvelmatu(i+1,k+1)/counter(i+1,k+1);
            outputvelmatv(i+1,k+1)=outputvelmatv(i+1,k+1)/counter(i+1,k+1);
            t=t+1;
        end
    end
end
clear counter;
outputvel(any(outputvel==0,2),:)=[];
outputvelmatu(any(outputvelmatu==0,2),:)=[];
outputvelmatv(any(outputvelmatv==0,2),:)=[];
toc

end

