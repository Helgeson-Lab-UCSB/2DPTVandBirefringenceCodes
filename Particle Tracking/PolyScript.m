%inputs---------------------------------------------
X=zeros(1,7);

%image parameters
lengthperpix=1/300;%distance per pixel
lengthunit='mm'; %unit of distance
framerate=280;%frames per time
timeunit='sec'; %unit of time
imagexres=1024;%pixel resolution of images
imageyres=544;
numbofim=500;%number of images
%imageloc='C:\Users\Patrick\12-17-18\Microfourrollmill Glycerol\Strain Rate';%image files location
imageloc='C:\Users\Patrick\Microfourrollmill Glycerol\Strain Rate';%image files location
imagename='0.2V Extensional';%image name without '0000','0001',etc or '.tif'

%graphing parameters
xcenter=331;%pixel center of geometry
ycenter=311;
maxMVGT=1000.0;%colorbar maximum for MVGT (1/timeunit)

%velocity and gradient quality parameters
trajectorygap=2; %images between velocity calculation
yres=1;%y slice in derivative calculation (please leave as 1)
xres=yres; %x slice in derivative calculation (please leave as 1)
p=30;
polyfitdeg=10;
roixmin=88;
roixmax=650;
roiymin=20;
roiymax=500;

%binning parameters
velbinsize=8;%velocity bin size (n x n square)
gradbinsize=8;%gradient bin size (n x n square)
edgecutoff=2;%bins on the edge cutoff

%PTV parameters
particlediameter=18;%maximum pixel diameter of tracer particle
minintensity=15;%minimum tracer particle intensity to be considered a feature
searchrad=18;%pixel radius to search for particle between frames
frameappear=5;%frames a particle must appear in to be considered a feature
framedisappear=5;%frames a particle can disappear for before forgotten about

%---------------------------------------------------
VelocityGradientTensorForScriptPoly
BeamVolAvgs

imagename='0.2V Extensional';%image name without '0000','0001',etc or '.tif'
VelocityGradientTensorForScriptPoly
BeamVolAvgs

imagename='0.3V Extensional';%image name without '0000','0001',etc or '.tif'
VelocityGradientTensorForScriptPoly
BeamVolAvgs

imagename='0.4V Extensional';%image name without '0000','0001',etc or '.tif'
VelocityGradientTensorForScriptPoly
BeamVolAvgs

imagename='0.5V Extensional';%image name without '0000','0001',etc or '.tif'
VelocityGradientTensorForScriptPoly
BeamVolAvgs

imagename='0.6V Extensional';%image name without '0000','0001',etc or '.tif'
VelocityGradientTensorForScriptPoly
BeamVolAvgs

imagename='0.8V Extensional';%image name without '0000','0001',etc or '.tif'
VelocityGradientTensorForScriptPoly
BeamVolAvgs

imagename='1.0V Extensional';%image name without '0000','0001',etc or '.tif'
VelocityGradientTensorForScriptPoly
BeamVolAvgs

imagename='1.2V Extensional';%image name without '0000','0001',etc or '.tif'
VelocityGradientTensorForScriptPoly
BeamVolAvgs

imagename='1.4V Extensional';%image name without '0000','0001',etc or '.tif'
VelocityGradientTensorForScriptPoly
BeamVolAvgs

imagename='1.6V Extensional';%image name without '0000','0001',etc or '.tif'
VelocityGradientTensorForScriptPoly
BeamVolAvgs

imagename='1.8V Extensional';%image name without '0000','0001',etc or '.tif'
VelocityGradientTensorForScriptPoly
BeamVolAvgs

imagename='2.0V Extensional';%image name without '0000','0001',etc or '.tif'
VelocityGradientTensorForScriptPoly
BeamVolAvgs
