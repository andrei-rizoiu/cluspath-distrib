function descriptiveGraphics( varargin )

	[filename, files, tags, vals, dtimestamp, clusassign, saveToFile, dumpEpsFiles] = process_options (varargin , ...
		'filename', [], 'files', [], 'tags', [], 'vals', [], 'dtimestamp', [], 'clusassign', [], 'saveToFile', false, ... 
		'dumpEpsFiles', false );
		
	if ( ~isempty(filename) )
		[files, tags, vals] = loadCsvFileMissing(filename, true);
		dtimestamp = vals(:, 1);
		vals(:, 1) = [];
		clusassign = vals(:, end);
		vals(:, end) = [];
	endif
	
	if ( isempty(files) || isempty(tags) || isempty(vals) || isempty(dtimestamp) || isempty(clusassign) )
		error("I don't have all the data I need to work!");
	endif
	
	noclusters = max(clusassign);	#number of clusters
	years = unique(dtimestamp)';	#vector of years
	entities = unique(files);		#an array with the list of entitites
	noObservations = length(files) / length(entities); #number of observations per entity
	noEntities = length(entities);	#number of entities in the dataset
	noYears = length(years);	#number of years in the dataset
	nfiles = zeros(1, length(files));	#nfiles will be the numeric representation of files over entities
	for i=1:noEntities
		pos = strcmp(files, entities{i});
		nfiles(pos) = i;
	endfor
	ystep = round(noEntities/10) -1;
	#remove underscores from entities names
	for i=1:noEntities
		entities{i} = strrep(entities{i}, '_', ' ');
	endfor
	colors={"blue", "cyan", "red", "magenta", "green", "yellow", "cyan", "blue", "black", "red", "magenta", "green", "yellow", "black", "blue", "cyan", "red", "magenta", "green", "yellow", "black", "blue", "cyan", "red", "magenta", "green", "yellow", "black"};
	markers={"+", "*", "o", "x", "^", "+", "*", "o", "*", "o", "+", "^", "*", "o","+", "*", "o", "x", "^", "+", "*", "o", "*", "o", "+", "^", "*", "o"};
	
	mainfig = figure(1);
	sz = size(clusassign);
	if (sz(2) == 1)	#column vector
		sz(2) = sz(1);
		sz(1) = 1;
	endif
	R = random("normal", 0, 0.05, sz );

	plots = [];
	legends = {};
	
	#################### START THE FIRST PLOT ####################################
	# first put the graphic with the observation distributions over time
	for i=1:noclusters
		x = zeros(1, length(dtimestamp));
		x(clusassign == i) = clusassign(clusassign == i);
		x(clusassign ~= i) = NaN;
		x = x + R;
		plots = [plots ; x];
		legends = [ legends {sprintf("\\mu_{%d}", i)}];
	endfor
	size(plots);
	subplot(2, 2, 1);
	h = plot(dtimestamp, plots, 'x');
	# set colors
	for i=1:noclusters
		set (h(i), "color", colors{i}, "marker", markers{i}, "linestyle", "none");
	endfor
	ylabel('Cluster number');
	xlabel('Time');
	title('Observations in clusters over time');
	if (noclusters <= 10)
		legend(legends, "location", "northwest");
	endif
	set(gca, 'Xtick', years(1:60:end));
	set(gca, 'XTickLabel', datestr(years(1:60:end)));
	set(gca, 'Ytick', [1:noclusters] );
	axis([min(years) max(years)+1 0.5 noclusters+0.5]);
	grid on;

	if (dumpEpsFiles)
		secfig = figure(mainfig+1);
		clf;
		
		h = plot(dtimestamp, plots, 'x');
		# set colors
		for i=1:noclusters
#			set (h(i), "color", colors{i}, "linewidth", 5, "marker", markers{i}, "markersize", 10);
			set (h(i), "color", colors{i}, "marker", markers{i}, "linestyle", "none");
		endfor
		set(gca, 'Xtick', years(1:60:end));
		set(gca, 'XTickLabel', datestr(years(1:60:end), "dd/mm/yy"));
		set(gca, 'Ytick', [1:noclusters] );
		axis([min(years) max(years)+1 0.5 noclusters+0.5]);
		ylabel('Cluster number', "fontsize", 20); #, "fontsize", 35
		xlabel('Time', "fontsize", 20); #, "fontsize", 35
		title('Observations in clusters over time', "fontsize", 20); #, "fontsize", 35
		legend(legends, "location", "northwest");
		legend("right");
		grid on;

		set (gca,'fontsize',15); # set (gca,'fontsize',30);

		print -tight -color -FHelvetica:15 -deps2 '01) graphic-observations-clusters-vs-time.eps' #-F:30
	endif
	
	################################# DONE WITH THE FIRST PLOT ######################################
	#go back to main figure
	figure(mainfig);
	#histograms over clusters over years
	subplot(2, 2, 2);
	bars = [];
	legends = {};
	for i=1:noclusters
		[x, foo] = hist(dtimestamp(clusassign == i), years);
		bars = [bars x'];
		legends = [ legends {sprintf("\\mu_{%d}", i)}];
	endfor
	h = bar(years, bars, 'stacked');
	# set colors
	for i=1:noclusters
		set (h(i), "facecolor", colors{i});
	endfor
	axis([min(years)-1 max(years)+1 0 noEntities]);
	set(gca, 'XTickLabel', datestr([years(1:60:end) max(years)+1], "dd/mm/yy") );
	title('Cluster distribution over years');
	xlabel('Years');
	ylabel('Cluster distribution');
	if (noclusters <= 10)
		legend(legends, "location", "eastoutside");
	endif
	
	if (dumpEpsFiles)
		secfig = figure(mainfig+1);
		clf;
		
		h = bar(years, bars, 'stacked');
		# set colors
		for i=1:noclusters
			set (h(i), "facecolor", colors{i});
		endfor
		axis([min(years)-1 max(years)+1]); #0 noEntities
		set(gca, 'XTickLabel', datestr([years(1:60:end) max(years)+1], "dd/mm/yy") );
		title('Cluster distribution over time', "fontsize", 20);
		xlabel('Time', "fontsize", 20);
		ylabel('Cluster distribution', "fontsize", 20);
		
		set (gca,'fontsize',15);
		print -tight -color -FHelvetica:15 -deps2 '02) graphic-clusters-vs-time.eps'	
	endif

	################################################# DONE WITH THE SECOND PLOT #########################################
	# let's skip it altogether	
	#go back to main figure
	figure(mainfig);
	#histograms over clusters over entities
	subplot(2, 2, 4);
	bars = [];
	legends = {};
	for i=1:noclusters
		[x, foo] = hist(nfiles(clusassign == i), [1:noEntities]);
		bars = [bars x'];
		legends = [ legends {sprintf("\\mu_{%d}", i)}];
	endfor
	h = barh([1:noEntities], bars, 'stacked');
	# set colors
	for i=1:noclusters
		set (h(i), "facecolor", colors{i});
	endfor
	axis([ 0 noObservations 0 noEntities+1]);
	title('Entity segmentation over time');
	xlabel('Time');
#	ylabel('Countries');

	set(gca, 'XTickLabel', datestr([years(1:60:end) max(years)+1], "dd/mm/yy") );
	
	set(gca,'Ytick',1:ystep:noEntities);
	set(gca, 'YTickLabel', entities(1:ystep:noEntities));
	if (noclusters <= 10)
		legend(legends, "location", "eastoutside");
	endif
	
	if (dumpEpsFiles)
		secfig = figure(mainfig+1);
		clf;
		
		h = barh([1:noEntities], bars, 'stacked');
		# set colors
		for i=1:noclusters
			set (h(i), "facecolor", colors{i});
		endfor
		axis([ 0 noObservations 0 noEntities+1]);
		title('Entity segmentation over time', "fontsize", 20);
		xlabel('Time', "fontsize", 20);

		set(gca, 'XTickLabel', datestr([years(1:60:end) max(years)+1], "dd/mm/yy") );
		set(gca,'Ytick',1:ystep:noEntities);
		set(gca, 'YTickLabel', entities(1:ystep:noEntities));
		
		set(gca,'fontsize',15);
		print -tight -color -FHelvetica:15 -deps2 '04) graphic-entity-segmentation-bars.eps'	
	endif
	
	###################################### DONE WITH THE THIRD PLOT #######################################"
	#go back to main figure
	figure(mainfig);
	# grafic for segmentation chunks for each country
	subplot(2, 2, 3);
	plots = [];
	for i=1:noclusters
		x = zeros(1, length(dtimestamp));
		x(clusassign == i) = nfiles(clusassign == i);
		x(clusassign ~= i) = NaN;
		plots = [plots ; x];
	endfor
	
	h = plot(dtimestamp, plots, '.');
	# set colors
	for i=1:noclusters
		set (h(i), "color", colors{i}, "linestyle", "none");
	endfor
	xlabel('Year');
#	ylabel('Entity segmentation');
	title('Segmentation for entities over time');
	
	axis([min(years)-1 max(years)+1 0 noEntities+1]);
	set(gca, 'Xtick', [min(years):10:max(years)+1] );
	set(gca,'Ytick',1:ystep:noEntities);
	set(gca, 'YTickLabel', entities(1:ystep:noEntities));
	if (noclusters <= 10)
		legend(legends, "location", "westoutside");
	endif
	
	if (dumpEpsFiles)
		secfig = figure(mainfig+1);
		clf;
		
		h = plot(dtimestamp, plots, '*');
		# set colors
		for i=1:noclusters
			set (h(i), "color", colors{i}, "linestyle", "none", "markersize", 2);
		endfor
		xlabel('Year', "fontsize", 20);
		title('Segmentation for entities over time', "fontsize", 20);
	
		axis([min(years)-1 max(years)+1 0 noEntities+1]);
		set(gca, 'Xtick', years(1:60:end));
		set(gca, 'XTickLabel', datestr(years(1:60:end), "dd/mm/yy"));
		set(gca,'Ytick',1:ystep:noEntities);
		set(gca, 'YTickLabel', entities(1:ystep:noEntities));
		
		set(gca,'fontsize',15);
		print -tight -FHelvetica:15 color -deps2 '03) graphic-entity-segmentation-points.eps'	
	endif
	
	#go back to main figure
	figure(mainfig);
	if (saveToFile)
		print('output-descriptive-graphics.png', '-dpng');
	endif

endfunction
