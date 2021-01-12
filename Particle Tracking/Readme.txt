This folder includes several codes for PTV and image analysis

Patrick Corona (corona@ucsb.edu)

-"BirefringenceCalc.m"
	-Calculates the redardance and direction of orientation from birefringence
	measurements (images) made through linear crossed polarizers at 0 and 45 degrees
	and the backgroiund

-"LOWESSScript.m/PolyScript.m"
	-Calculates the 2D velocity gradient tensor from a series of images from movies
	-If too slow, PolyScript uses a high order polynominal fit rather than LOWESS
	which is significantly faster and provides similar results for the conditions I tested
	-Script utilizes many functions from the Nfeature3 folder, so ensure this is 
	added to the path before running

-"PlotGrad.m", "QuickRebin.m", "ElasticityPaperCalc.m"
	-Some scripts I used for quickly calculating new results from LOWESSScript and PolyScript

-Nfeature3
	-Contains many scripts and functions that are called in the above codes