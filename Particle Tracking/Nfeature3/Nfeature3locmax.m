function [loc x y] = Nfeature3locmax(img,lambda,w,field,Imin)

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

