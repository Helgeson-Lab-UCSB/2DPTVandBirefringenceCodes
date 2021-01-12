%inputs---------------------------------------------
X=zeros(1,7);

%image parameters
lengthperpix=2/187;%distance per pixel
lengthunit='mm'; %unit of distance
framerate=280;%frames per time
timeunit='sec'; %unit of time
imagexres=1024;%pixel resolution of images
imageyres=544;
numbofim=200;%number of images
%imageloc='C:\Users\Patrick\12-17-18\Microfourrollmill Glycerol\Strain Rate';%image files location
imageloc='C:\Users\Patrick\19-1-29\fformsaxs 0.1wt% cnf';%image files location
imagename='fformsaxs 0.1wt%cnf q2=0.1 frr=0 280fps(2)';%image name without '0000','0001',etc or '.tif'

%graphing parameters
xcenter=611;%pixel center of geometry
ycenter=254;
maxMVGT=100.0;%colorbar maximum for MVGT (1/timeunit)

%velocity and gradient quality parameters
trajectorygap=2; %images between velocity calculation
yres=1;%y slice in derivative calculation (please leave as 1)
xres=yres; %x slice in derivative calculation (please leave as 1)
p=50;
polyfitdeg=10;
roixmin=200;
roixmax=850;
roiymin=20;
roiymax=500;


%binning parameters
velbinsize=8;%velocity bin size (n x n square)
gradbinsize=8;%gradient bin size (n x n square)
edgecutoff=1;%bins on the edge cutoff

%PTV parameters
particlediameter=10;%maximum pixel diameter of tracer particle
minintensity=30;%minimum tracer particle intensity to be considered a feature
searchrad=4;%pixel radius to search for particle between frames
frameappear=5;%frames a particle must appear in to be considered a feature
framedisappear=10;%frames a particle can disappear for before forgotten about

%---------------------------------------------------
VelocityGradientTensorForScript
BeamVolAvgs

imagename='fformsaxs 0.1wt%cnf q2=0.1 frr=-1 280fps(3)';%image name without '0000','0001',etc or '.tif'
xcenter=494;%pixel center of geometry
ycenter=253;
VelocityGradientTensorForScript
BeamVolAvgs

imagename='fformsaxs 0.1wt%cnf q2=0.1 frr=1 280fps(3)';%image name without '0000','0001',etc or '.tif'
xcenter=494;%pixel center of geometry
ycenter=253;
VelocityGradientTensorForScript
BeamVolAvgs
