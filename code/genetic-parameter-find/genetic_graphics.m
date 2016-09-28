function genetic_graphics(History, bestIndividual)

	if ( isempty(History))
		error("I cannot work with an empty History. Give me the data structure or the file in which it is saved.");
	endif

	if (nargin < 2)
		bestIndividual = [];
	endif

	# if here, we were given a History to continue
	# is it a structure (an Octave object?)
	if (isstruct(History))
		ng = length(History);
		Configs = History(ng).Configs;
		Outputs = History(ng).Outputs;
		Measures = History(ng).Measures;
	else
		# maybe it is a file to load from; lets see what is in it?
		[ vars ] = who("-file", History);
		correct = length(vars) == 1;	%need to have the 1 variables here after
		for i=1:length(vars)
			if ( strcmp( vars{i}, "History" ) == 0 )
				correct = false;
			endif
		endfor
	
		if ( correct )
			load( History );
			ng = length(History);
			Configs = History(ng).Configs;
			Outputs = History(ng).Outputs;
			Measures = History(ng).Measures;
		else
			error(sprintf("The given file '%s' is not a valid history file!"));
		endif
	endif
	
	# plot the MDvar-Tvar paretor front
	close all;
	mainfig = figure(1);
	subplot(2, 2, 1);
	legends = plot_graphic_2measures(History, 1, 2, bestIndividual, false);
	
	# plot the MDvar-ShaP pareto front
	figure(mainfig);
	subplot(2, 2, 2);
	plot_graphic_2measures(History, 1, 3, bestIndividual, false);
		
	# plot the ShaP-SPass pareto front
	figure(mainfig);
	subplot(2, 2, 3);
	plot_graphic_2measures(History, 2, 3, bestIndividual, false);
	
	# plot the MDvar-SPass pareto front
	figure(mainfig);
	subplot(2, 2, 4);
	plot_graphic_2measures(History, 1, 4, bestIndividual, false);
	
	## put the legend outside
	## dummy axes
	## styles to be used
	tmp = axes (); hold on;
	dumx = -1000;
	dumy = -1000;

	# plot dummy data for legend
	plot(dumx, dumy, "^r", 'LineWidth', 3, 'markersize', 3);
	plot(dumx, dumy, '*b', 'LineWidth', 4, 'markersize', 5); 
	plot( [dumx-10 ; dumx], [dumy-10 ; dumy], "-b");
	
	legend ( {"Final generation", "Chosen solution", "Pareto Front"}, ... 
		"location", "southoutside", "orientation", "horizontal");
	legend("boxon"); 
	axis ([-1 100 -10 2])
	axis ("off")  ## hides the axes lines only

	
	print("genetic-evolution-graphics.pdf", "-dpdf", "-FHelvetica:11", "-tight");
	system("pdfcrop genetic-evolution-graphics.pdf ; mv genetic-evolution-graphics-crop.pdf genetic-evolution-graphics.pdf");
endfunction

