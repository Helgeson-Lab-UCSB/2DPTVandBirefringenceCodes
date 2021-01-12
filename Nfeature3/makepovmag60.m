function d = makepovmag60(pts,filename)
% makepov.m		5-12-03    Maria Kilfoil

% ; pts contains the 3D points
% ; zmag stretches in z-direction   (example:  zmag=2)
% ; radius changes size of particles (in same units as original data)
% ; /bw for black & white
% ; /nobox to remove outer box
% ; margin = 'whatever' to add a little margin to the box (make the box bigger)
% ; camera = [x,y,z] to relocate camera
% ; lookat = [x,y,z] to relocate where camera looks at
% ; color = [r,g,b] data for each point

% 
zmag=1;
% if (not keyword_set(bw)) then bw=0 else bw=1
% if (not keyword_set(radius)) then radius=0.5
% 
tr0=pts;
% 
tr0(:,3)=tr0(:,3)+20.0;
tr0(:,3)=tr0(:,3)*zmag;
% avgstd,tr0(0,*),/quiet,xres
% avgstd,tr0(1,*),/quiet,yres
% avgstd,tr0(2,*),/quiet,zres
% zres(2,*)=zres(2,*) > 0
% dz=zres(3,*)-zres(2,*)
% zres(2,*)=zres(3,*)+dz*zmag
% zres(0,*)=zres(3,*)+dz*0.5
% 
xres(2)=-0.6418340;
yres(2)=-0.483080;
zres(2)=68.603;
xres(3)=51.2022;
% yres(3)=70.1979;
yres(3)=35.6279;
zres(3)=19.8447;
% 
% close,1
% openw,1,filename
fid = fopen(filename,'w');

fprintf(fid,'#include "colors.inc"\n');
fprintf(fid,'#include "shapes.inc"\n');
fprintf(fid,'#include "textures.inc"\n');
% xbar=(xres(3)-xres(2))*0.5+xres(2)
% ybar=(yres(3)-yres(2))*0.5+yres(2)
% zbar=(zres(3)-zres(2))*0.5+zres(2)
% ; print,xbar,ybar,zbar

fprintf(fid,'camera {\n');
fprintf(fid,'   location <      165.5870      45.0405      1270.118>\n');
fprintf(fid,'   look_at <  55.5870   45.0405  25.7236>\n');
fprintf(fid,'   angle 6.0 }\n');

fprintf(fid,'light_source { <25,52,135> color White}\n');
fprintf(fid,'background {color rgb <1,1,1>}\n');
fprintf(fid,'#declare R1 = 0.47;\n');
fprintf(fid,'#declare R2 = 0.12;\n');
fprintf(fid,'#declare ballfin = finish {\n');
fprintf(fid,'   reflection 0.0 diffuse 0.4 ambient 0.7 phong 1.0 phong_size 300}\n');

fprintf(fid,'#declare redcolor = texture {\n');
fprintf(fid,'   pigment {color rgb <1,0,0>}\n');
fprintf(fid,'   finish { ballfin }}\n');
fprintf(fid,'#declare boxtex = texture {\n');
fprintf(fid,'   pigment {color rgb <0.6,0.8,1>}\n');
fprintf(fid,'   finish { reflection 0.0\n');
fprintf(fid,'      diffuse 0.4 ambient 0.7 phong 0.0 }}\n');

nw1=length(tr0(:,1));
for i=1:nw1,
% 	if (not keyword_set(color)) then begin
		redsphere(fid,tr0(i,1:3));
% 	endif else begin
% 		colsphere,tr0(*,i),reform(color(*,i))
% 	endelse
end
% ; ============================================================

% if (keyword_set(margin)) then begin
% 	xres(2)=xres(2)-margin
% 	yres(2)=yres(2)-margin
% 	zres(2)=zres(2)-margin
% 	xres(3)=xres(3)+margin
% 	yres(3)=yres(3)+margin
% 	zres(3)=zres(3)+margin
% endif

	p0=[xres(2),yres(2),zres(2)]; p1=[xres(3),yres(2),zres(2)];
	colcylinder(fid,p0,p1);
    p0=[xres(2),yres(3),zres(2)]; p1=[xres(3),yres(3),zres(2)];
	colcylinder(fid,p0,p1);
    p0=[xres(2),yres(3),zres(3)]; p1=[xres(3),yres(3),zres(3)];
	colcylinder(fid,p0,p1);
	p0=[xres(2),yres(2),zres(3)]; p1=[xres(3),yres(2),zres(3)];
	colcylinder(fid,p0,p1);
	p0=[xres(2),yres(2),zres(2)]; p1=[xres(2),yres(3),zres(2)];
	colcylinder(fid,p0,p1);
	p0=[xres(3),yres(2),zres(2)]; p1=[xres(3),yres(3),zres(2)];
	colcylinder(fid,p0,p1);
	p0=[xres(3),yres(2),zres(3)]; p1=[xres(3),yres(3),zres(3)];
	colcylinder(fid,p0,p1);
	p0=[xres(2),yres(2),zres(3)]; p1=[xres(2),yres(3),zres(3)];
	colcylinder(fid,p0,p1);
	p0=[xres(2),yres(2),zres(2)]; p1=[xres(2),yres(2),zres(3)];
	colcylinder(fid,p0,p1);
	p0=[xres(3),yres(2),zres(2)]; p1=[xres(3),yres(2),zres(3)];
	colcylinder(fid,p0,p1);
	p0=[xres(3),yres(3),zres(2)]; p1=[xres(3),yres(3),zres(3)];
	colcylinder(fid,p0,p1);
	p0=[xres(2),yres(3),zres(2)]; p1=[xres(2),yres(3),zres(3)];
	colcylinder(fid,p0,p1);

st = fclose(fid)
