rr=[]
for i=1:60
    im=imread(strcat('C:\Users\HelgesonGp\Desktop\Calibration with Karo Light Corn Syrup\Overlap of Multiple Measurement\Measurement1\Measurement10',num2str(i-1,'%03.0f'),'.tif'));
    r=Nfeature3(im,1,3,2,40);
    temp=r;
    temp(:,6)=i;
    rr=[rr;temp];
end

ss=rr(:,1:2);
ss(:,3)=rr(:,6);

[lub]=trackmem(ss,5,2,5,5);

j=max(lub(:,4));
for k=1:j
    pt=lub(find(lub(:,4)==k),:);
    pt=sortrows(pt,3);
    hold on;
    plot(pt(:,1),pt(:,2),'b')
end
box on;

start=lub(find(lub(:,3)==1),:);
finish=lub(find(lub(:,3)==15),:);
[id ia ib]=intersect(start(:,4),finish(:,4),'rows');
pt1=start(ia,:);
pt2=finish(ib,:);
a=size(pt1);

figure
for m=1:a(1)
    vectarrow(pt1(m,1:2),pt2(m,1:2));
    hold on;
end