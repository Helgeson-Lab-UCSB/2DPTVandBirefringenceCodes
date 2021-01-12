a=size(im);
avg=[];
for i=1:a(2)
    avg(i,1)=i;
    avg(i,2)=mean(im(:,i));
end