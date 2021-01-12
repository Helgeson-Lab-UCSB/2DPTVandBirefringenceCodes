function [r] = Nfeature3rawimtest(img,lambda,w, Imin)
%      7-29-03  Maria Kilfoil
%extent should be 2*w+1 in which w is the same as bpass2.
% 	Finds and measures roughly circular 'features' within an image.
%  CALLING SEQUENCE:
% 	f = feature( image, diameter,field,Imin )
%  INPUTS:
% 	image:	(nx,ny) array which presumably contains some features worth finding
% 	extent: a parameter which should be a little greater than the diameter of the 
%           largest features in the image. Diameter MUST BE ODD valued.
% 	separation: an optional parameter which specifies the minimum allowable separation 
%           between feature centers. The default value is diameter+1.
% 	masscut: (*to add*) Setting this parameter saves runtime by reducing the runtime wasted on 
%           low mass 'noise' features.
% 	Imin: 	Set this optional parameter to the minimum allowed value for the peak 
%           brightness of a feature. Useful for limiting the number of spurious features in
% 	noisy images.
% 	field: 	Set this keyword if image is actually just one field of an interlaced 
%           (e.g. video) image. All the masks will then be constructed with a 2:1 aspect ratio.
% 	iterate: (to add later) if the refined centroid position is too far from the initial estimate, 
%           iteratively recalc. the centroid using the last cetroid to position the mask.  This 
% 	can be useful for really noisy data, or data with flat (e.g. saturated) peaks.  
%           Use with caution- it may 'climb' hills and give you multiple hits.
%  OUTPUTS:
% 		f(:,1):	the x centroid positions, in pixels.
% 		f(:,2): the y centroid positions, in pixels. 
% 		f(:,3): integrated brightness of the features. ("mass")
% 		f(:,4): the square of the radius of gyration of the features.
% 		    (second moment of the "mass" distribution, where mass=intensity)
% 		f(:,5): eccentricity, which should be zero for circularly symmetric features and 
%                   order one for very elongated images.
%  RESTRICTIONS:
%       To work properly, the image must consist of bright, circularly symmetric regions 
%       on a roughly zero-valued background. To find dark features, the image should be 
%       inverted and the background subtracted. If the image contains a large amount of 
%       high spatial frequency noise, performance will be improved by first filtering the image.
%       BPASS will remove high spatial frequency noise, and subtract the image background. 
%       Individual features should NOT overlap.
%
extent=2*w+1;
image = bpass3(img,lambda,w);
% image=img;
if (mod(extent,2) == 0),
    disp('Requires an odd extent.  Adding 1...');
    extent = extent + 1;
end
sz = size(image);
nx = sz(2);
ny = sz(1);
% if n_params() eq 2 then sep = extent+1
sep = extent; %sep separation   


%       Put a border around the image to prevent mask out-of-bounds
% a = zeros( ny + extent, nx + extent );
% a(fix(extent/2)+1:(fix(extent/2))+ny,fix(extent/2)+1:(fix(extent/2))+nx) = image;
%New   It seems that it has already had a border of w width
a=image;
whos a
% nx = nx + extent;   % a doesn't change at all if we don't add a border to
% it


%       Finding local maxima
field=2;
% loc= localmax3nodilate(image,sep,field,Imin);
loc= localmax3(image,sep,field,Imin);
if (loc(1) == -1)
    r = -1;
    return
else
    ;
end
y = mod(loc,ny);
x = fix(loc / ny+1);  %the pixel position of the maximum
xraw=x;
yraw=y;
% ** make x,y double precision? **
nmax=length(loc);           %number of feathers?
disp(strcat(num2str(nmax,'%01.0f'),' features found.'));
% image = bpass3(img,lambda,w);
% a=image;
% img=double(img);
% img=img/(max(max(img))-min(min(img)))*255;
% image=double(img);

%	Setup some result arrays
m=zeros(nmax,1);
xc = zeros(nmax,1);
yc = zeros(nmax,1);
rg = zeros(nmax,1);
e  = zeros(nmax,1);
xn = zeros(nmax,1);
yn = zeros(nmax,1);
masksz=zeros(nmax,1);
%	Estimate the mass	
for k=1:nmax,
    edge_search=1;
    wmask=0;
     while edge_search==1
      wmask=wmask+1;        
      extent=2*wmask+1;
      xl = x(k) - fix(extent/2);
      xh = xl + extent -1;         %xh-xl=extent-1 which is diameter of the feature
      yl = y(k) - fix(extent/2);
      yh = yl + extent -1;           %one diameter high
      subimg=double(a(fix(yl):fix(yh),fix(xl):fix(xh)));
%       reach_edge=sum(a(y(k),fix(xl):fix(xh))<10)+ sum(a(fix(yl):fix(yh),x(k))<10);
      if sum(sum(subimg<10))>0
          edge_search=0;
      end
     end
%      extent=extent+2;
      xl = x(k) - fix(extent/2);
      xh = xl + extent -1;         %xh-xl=extent-1 which is diameter of the feature
      yl = y(k) - fix(extent/2);
      yh = yl + extent -1;           %one diameter high
     subimg=double(image(fix(yl):fix(yh),fix(xl):fix(xh)));
     masksz(k)=extent;
%       Set up some masks
    rsq = rsqd( extent,extent );
    t = thetarr( extent );
mask = le(rsq,(extent/2)^2);       %get a circular unit mask
mask2 = ones(1,extent)'*[1:extent]; 
mask2 = mask2.*mask;              %get a circular mask, increasing with collum
mask3= (rsq.*mask) + (1/6);
cen = (extent-1)/2 +1;                %mask center
% cmask = vpa(cos(sym('2')*t)).*mask;  
% ultra high presision, since Matlab and IDL differ here
cmask = cos(2*t).*mask;
% smask = vpa(sin(sym('2')*t)).*mask;
smask = sin(2*t).*mask;
cmask(cen,cen) = 0.0;
smask(cen,cen) = 0.0;   %give zero value to the center of the mask
    suba = zeros(extent, extent, nmax);
    xmask = mask2;
    ymask = mask2';
	ycen = cen;                   %just (extent+1)/2
    %get total mass for each feature
    m(k) = sum(sum(subimg.*mask));   
    %	Calculate feature centers
    xc(k) = sum(sum(subimg.*xmask));  
	yc(k) = sum(sum(subimg.*ymask));
    %	Correct for the 'offset' of the centroid masks
    xc(k) = xc(k)/m(k) - (extent+1)/2;             
    yc(k) = yc(k)/m(k) - (extent+1)/2;  %get mass center
    %	Update the positions and correct for the width of the 'border'
    if abs(xc(k))>=1 | abs(yc(k))>=1
    x(k)=x(k)+fix(xc(k));
    xc(k)=xc(k)-fix(xc(k));
    y(k)=y(k)+fix(yc(k));
    yc(k)=yc(k)-fix(xc(k));
    xl = x(k) - fix(extent/2);
    xh = xl + extent -1;         %xh-xl=extent-1 which is diameter of the feature
    yl = y(k) - fix(extent/2);
    yh = yl + extent -1;           %one diameter high
    end
%     x(k) = x(k) + xc(k) - 0*fix(extent/2);
%     y(k) = y(k) + yc(k) - 0*fix(extent/2);
    xn(k)=x(k);
    yn(k)=y(k);
%%%%%%Fracshift; refinement to obtain subpixel resolution%%%%%%%%%%%%%%%%%%%
    xcn(k)=xc(k);
    ycn(k)=yc(k);
%	Construct the subarray and calculate the mass, radii of gyration squared, eccentricity
    nshift=1;
    while nshift==1
       xn(k)=x(k)+xcn(k);
       yn(k)=y(k)+ycn(k);
       if fix(yl)-1>0 & fix(yh)+1<=ny & fix(xl)-1>0 & fix(xh)+1<=nx
    suba = fracshift( double(image(fix(yl)-1:fix(yh)+1,fix(xl)-1:fix(xh)+1)), -xcn(k) , -ycn(k) );
    subaa=suba(2:end-1, 2:end-1);
       else
    suba = fracshift( double(image(fix(yl):fix(yh),fix(xl):fix(xh))), -xcn(k) , -ycn(k) );  
    subaa=suba;
       end
    m(k) = sum(sum(( subaa.*mask )));             % mass
    rg(k) = (sum(sum( subaa.*mask3 ))) / m(k);    % radii of gyration squared
    tmp = sqrt(( (sum(sum( subaa.*cmask )))^2 ) +( (sum(sum( subaa.*smask )))^2 )); 
    tmp2 = (m(k)-subaa(cen,ycen)+1e-6);
    e(k) = tmp/tmp2;    % eccentricity
    xc(k) = sum(sum(double(subaa).*xmask));  
	yc(k) = sum(sum(double(subaa).*ymask));
    xc(k) = xc(k)/ m(k) - ((double(extent(1))+1.0)/2.0);
    yc(k) = yc(k)/ m(k) - ((double(extent(1))+1.0)/2.0);
    xcn(k)=xc(k)+xcn(k);
    ycn(k)=yc(k)+ycn(k);
    if abs(xcn(k))>=1 | abs(ycn(k))>=1
    x(k)=x(k)+fix(xcn(k));
    xcn(k)=xcn(k)-fix(xcn(k));
    y(k)=y(k)+fix(ycn(k));
    ycn(k)=ycn(k)-fix(xcn(k));
    xl = x(k) - fix(extent/2);
    xh = xl + extent -1;         %xh-xl=extent-1 which is diameter of the feature
    yl = y(k) - fix(extent/2);
    yh = yl + extent -1;           %one diameter high
    end
%     xn(k) = x(k);
%     yn(k) = y(k);
%     x(k) = x(k) + xc(k) ;
%     y(k) = y(k) + yc(k) ;
    if xc(k)<0.0001 & yc(k)<0.0001
        nshift=0;
    end
    end 
end

% %%%%%%%%%%%shrink mask size for features shifted more than a pixel%%%
% id=find(abs(xcn)>1 | abs(ycn)>1); whos id
% x(id)=xraw(id);
% y(id)=yraw(id);
% for k=1:length(id),
%      extent=masksz(id(k))-4; 
%      if extent<1
%          extent=1;
%      end
%      masksz(id(k))=extent;
%       xl = x(id(k)) - fix(extent/2);
%       xh = xl + extent -1;         %xh-xl=extent-1 which is diameter of the feature
%       yl = y(id(k)) - fix(extent/2);
%       yh = yl + extent -1;           %one diameter high
%      subimg=double(a(fix(yl):fix(yh),fix(xl):fix(xh)));
% %       Set up some masks
%     rsq = rsqd( extent,extent );
%     t = thetarr( extent );
% mask = le(rsq,(extent/2)^2);       %get a circular unit mask
% mask2 = ones(1,extent)'*[1:extent]; 
% mask2 = mask2.*mask;              %get a circular mask, increasing with collum
% mask3= (rsq.*mask) + (1/6);
% cen = (extent-1)/2 +1;                %mask center
% % cmask = vpa(cos(sym('2')*t)).*mask;  
% % ultra high presision, since Matlab and IDL differ here
% cmask = cos(2*t).*mask;
% % smask = vpa(sin(sym('2')*t)).*mask;
% smask = sin(2*t).*mask;
% cmask(cen,cen) = 0.0;
% smask(cen,cen) = 0.0;   %give zero value to the center of the mask
%     suba = zeros(extent, extent, nmax);
%     xmask = mask2;
%     ymask = mask2';
% 	ycen = cen;                   %just (extent+1)/2
%     %get total mass for each feature
%     m(id(k)) = sum(sum(subimg.*mask));   
%     %	Calculate feature centers
%     xc(id(k)) = sum(sum(subimg.*xmask));  
% 	yc(id(k)) = sum(sum(subimg.*ymask));
%     %	Correct for the 'offset' of the centroid masks
%     xc(id(k)) = xc(id(k))/m(id(k)) - (extent+1)/2;             
%     yc(id(k)) = yc(id(k))/m(id(k)) - (extent+1)/2;  %get mass center
%     %	Update the positions and correct for the width of the 'border'
%     x(id(k)) = x(id(k)) + xc(id(k)) - 0*fix(extent/2);
%     y(id(k)) = y(id(k)) + yc(id(k)) - 0*fix(extent/2);
%     xn(id(k))=x(id(k));
%     yn(id(k))=y(id(k));
% %%%%%%Fracshift; refinement to obtain subpixel resolution%%%%%%%%%%%%%%%%%%%
%     xcn(id(k))=xc(id(k));
%     ycn(id(k))=yc(id(k));
% %	Construct the subarray and calculate the mass, radii of gyration squared, eccentricity
%     nshift=1;
%     while nshift==1
%        if fix(yl)-1>0 & fix(yh)+1<=ny & fix(xl)-1>0 & fix(xh)+1<=nx
%     suba = fracshift( double(image(fix(yl)-1:fix(yh)+1,fix(xl)-1:fix(xh)+1)), -xcn(id(k)) , -ycn(id(k)) );
%     subaa=suba(2:end-1, 2:end-1);
%        else
%     suba = fracshift( double(image(fix(yl):fix(yh),fix(xl):fix(xh))), -xcn(id(k)) , -ycn(id(k)) );  
%     subaa=suba;
%        end
%     m(id(k)) = sum(sum(( subaa.*mask )));             % mass
%     rg(id(k)) = (sum(sum( subaa.*mask3 ))) / m(id(k));    % radii of gyration squared
%     tmp = sqrt(( (sum(sum( subaa.*cmask )))^2 ) +( (sum(sum( subaa.*smask )))^2 )); 
%     tmp2 = (m(id(k))-subaa(cen,ycen)+1e-6);
%     e(id(k)) = tmp/tmp2;    % eccentricity
%     xc(id(k)) = sum(sum(double(subaa).*xmask));  
% 	yc(id(k)) = sum(sum(double(subaa).*ymask));
%     xc(id(k)) = xc(id(k))/ m(id(k)) - ((double(extent(1))+1.0)/2.0);
%     yc(id(k)) = yc(id(k))/ m(id(k)) - ((double(extent(1))+1.0)/2.0);
%     xcn(id(k))=xc(id(k))+xcn(id(k));
%     ycn(id(k))=yc(id(k))+ycn(id(k)); 
%     xn(id(k)) = x(id(k));
%     yn(id(k)) = y(id(k));
%     x(id(k)) = x(id(k)) + xc(id(k)) ;
%     y(id(k)) = y(id(k)) + yc(id(k)) ;
%     if xc(id(k))<0.0001 & yc(id(k))<0.0001
%         nshift=0;
%     end
%     end 
% end

% id=find(abs(xcn)>1 | abs(ycn)>1); whos id
% x(id)=xraw(id);
% y(id)=yraw(id);
% for k=1:length(id),
%      extent=masksz(id(k))-2; 
%      if extent<1
%          extent=1;
%      end
%      masksz(id(k))=extent;
%       xl = x(id(k)) - fix(extent/2);
%       xh = xl + extent -1;         %xh-xl=extent-1 which is diameter of the feature
%       yl = y(id(k)) - fix(extent/2);
%       yh = yl + extent -1;           %one diameter high
%      subimg=double(a(fix(yl):fix(yh),fix(xl):fix(xh)));
% %       Set up some masks
%     rsq = rsqd( extent,extent );
%     t = thetarr( extent );
% mask = le(rsq,(extent/2)^2);       %get a circular unit mask
% mask2 = ones(1,extent)'*[1:extent]; 
% mask2 = mask2.*mask;              %get a circular mask, increasing with collum
% mask3= (rsq.*mask) + (1/6);
% cen = (extent-1)/2 +1;                %mask center
% % cmask = vpa(cos(sym('2')*t)).*mask;  
% % ultra high presision, since Matlab and IDL differ here
% cmask = cos(2*t).*mask;
% % smask = vpa(sin(sym('2')*t)).*mask;
% smask = sin(2*t).*mask;
% cmask(cen,cen) = 0.0;
% smask(cen,cen) = 0.0;   %give zero value to the center of the mask
%     suba = zeros(extent, extent, nmax);
%     xmask = mask2;
%     ymask = mask2';
% 	ycen = cen;                   %just (extent+1)/2
%     %get total mass for each feature
%     m(id(k)) = sum(sum(subimg.*mask));   
%     %	Calculate feature centers
%     xc(id(k)) = sum(sum(subimg.*xmask));  
% 	yc(id(k)) = sum(sum(subimg.*ymask));
%     %	Correct for the 'offset' of the centroid masks
%     xc(id(k)) = xc(id(k))/m(id(k)) - (extent+1)/2;             
%     yc(id(k)) = yc(id(k))/m(id(k)) - (extent+1)/2;  %get mass center
%     %	Update the positions and correct for the width of the 'border'
%     x(id(k)) = x(id(k)) + xc(id(k)) - 0*fix(extent/2);
%     y(id(k)) = y(id(k)) + yc(id(k)) - 0*fix(extent/2);
%     xn(id(k))=x(id(k));
%     yn(id(k))=y(id(k));
% %%%%%%Fracshift; refinement to obtain subpixel resolution%%%%%%%%%%%%%%%%%%%
%     xcn(id(k))=xc(id(k));
%     ycn(id(k))=yc(id(k));
% %	Construct the subarray and calculate the mass, radii of gyration squared, eccentricity
%     nshift=1;
%     while nshift==1
%        if fix(yl)-1>0 & fix(yh)+1<=ny & fix(xl)-1>0 & fix(xh)+1<=nx
%     suba = fracshift( double(image(fix(yl)-1:fix(yh)+1,fix(xl)-1:fix(xh)+1)), -xcn(id(k)) , -ycn(id(k)) );
%     subaa=suba(2:end-1, 2:end-1);
%        else
%     suba = fracshift( double(image(fix(yl):fix(yh),fix(xl):fix(xh))), -xcn(id(k)) , -ycn(id(k)) );  
%     subaa=suba;
%        end
%     m(id(k)) = sum(sum(( subaa.*mask )));             % mass
%     rg(id(k)) = (sum(sum( subaa.*mask3 ))) / m(id(k));    % radii of gyration squared
%     tmp = sqrt(( (sum(sum( subaa.*cmask )))^2 ) +( (sum(sum( subaa.*smask )))^2 )); 
%     tmp2 = (m(id(k))-subaa(cen,ycen)+1e-6);
%     e(id(k)) = tmp/tmp2;    % eccentricity
%     xc(id(k)) = sum(sum(double(subaa).*xmask));  
% 	yc(id(k)) = sum(sum(double(subaa).*ymask));
%     xc(id(k)) = xc(id(k))/ m(id(k)) - ((double(extent(1))+1.0)/2.0);
%     yc(id(k)) = yc(id(k))/ m(id(k)) - ((double(extent(1))+1.0)/2.0);
%     xcn(id(k))=xc(id(k))+xcn(id(k));
%     ycn(id(k))=yc(id(k))+ycn(id(k)); 
%     xn(id(k)) = x(id(k));
%     yn(id(k)) = y(id(k));
%     x(id(k)) = x(id(k)) + xc(id(k)) ;
%     y(id(k)) = y(id(k)) + yc(id(k)) ;
%     if xc(id(k))<0.0001 & yc(id(k))<0.0001
%         nshift=0;
%     end
%     end 
% end
% id=find(abs(xcn)>1 | abs(ycn)>1); whos id

r = [xn,yn,m,rg,e,masksz];
fshift=[xcn' ycn'];
%b;

