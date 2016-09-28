## Copyright (C) 2009 - 2013 Marian-Andrei Rizoiu.
##
## This file is part of the prototype I'm creating in
## preparation of my PhD thesis with the 
## Laboratory ERIC, University Lyon 2.
##
## If at any point, somebody finds this useful, he if free
## to use it. It is a free software; you can redistribute it 
## and/or modify it under the terms of the GNU General Public
## License as published by the Free Software Foundation;
## either version 3 of the License, or (at your option) any
## later version <http://www.gnu.org/licenses/>.

## Usage: temporalAwareMeasureColorMap('Alpha', 0, 'printToScreen', true, 'dumpEpsFiles', true )
##
## Creates the colored surface graphic for the temporal-aware dissimilarity matrix for different values of alpha

## Created: 19/02/2013 Marian-Andrei RIZOIU
##
function temporalAwareMeasureColorMap( varargin )

	[Alpha, printToScreen, dumpEpsFiles] = process_options (varargin , 'Alpha', 0, 'printToScreen', true, 'dumpEpsFiles', true );

	# calculate gammaD and gammaT
	gammaD = 1 + Alpha ; 
	if (gammaD > 1)
		gammaD = 1;
	endif

	gammaT = 1 - Alpha ; 
	if (gammaT > 1)
		gammaT = 1;
	endif

	# create the dx and dt vectors between 0 and 1 (already scaled) and calculate the values of the measure
	dx=dt=[0:0.005:1];
	[xx, yy] = meshgrid(dx, dt);
	zz = 1 - ( 1 - gammaD*(xx.^2)).*(1 - gammaT*(yy.^2));

	# lets create the graphic
	h=figure(1);
	clf;
	surface(dx, dt, zz); 
	shading flat;

	FN = findall(h,'-property','FontName');
	set(FN,'FontName', 'LiberationSans');
	FS = findall(h,'-property','FontSize');
	set(FS,'FontSize', 12);

	set(gca, 'Xtick', [0:0.1:1] );
	set(gca, 'Ytick', [0:0.1:1] );

	box("off");
	colorbar("EastOutside");
	
	title("Temporal-Aware Dissimilarity Measure ", "fontsize", 16); #, "fontsize", 40
	xlabel("Multidimensional distance", "fontsize", 15); #, "fontsize", 40
	ylabel("Temporal distance", "fontsize", 15); #, "fontsize", 40

	if (dumpEpsFiles)
		W = 7.08661417; # 18cm
		H = 4,72440945; # 12cm
		set(h,'PaperOrientation','portrait');
		set(h,'PaperSize',[H,W]);
		set(h,'PaperPosition',[0,0,W,H]);

		print -tight -color -deps2 -FLiberationSans 'surface_temporal_aware_measure.eps'	%save as eps and color -FHelvetica:14 
	endif
endfunction
