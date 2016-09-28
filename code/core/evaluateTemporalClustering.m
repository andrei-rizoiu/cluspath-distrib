## Usage: evaluateTemporalClustering('inputfile', 'temporal-average.csv', 'clusassign', clusassign, 'centroids',centroids, 'ctimestamp', ctimestamp);
##
## Given the input filename and a cluster assignement, this function calculates different 
## evaluation measures and prints them to the outputfid (separated by TAB). If 'printHeader' is true
## then only the header with measure names is printed. 'outputfid' is the stream fid to which to write.
##
## Parameters:
##	'inputfile' - the input dataset the clustering was performed on to work with; 
##	'clusassign' - the cluster assignment of the observations to clusters;
##	'centroids' - the centroids resulted from the clustering;
##	'ctimestamp' - the centroids timestamp;
##	'printHeader' (optional, default false) - if true, only prints header with the measures if calculates and exits;
##	'outputfid' (optional, default 1) - unix file id to write to. Default is 1 (on screen);
##	'algoname' (optional, default '') - prefix line of measures with this (usefull to output parameters etc.);
##	'algoheader' (optional, default '') - prefix for when printing only headers;
##	'outputfile' (optional, default []) - append to this file instead of outputfid. Cannot be used together with outputfid;
##	'normalize' (optional, default true) - if true, normalise input data before work;
##	'stardardize' (optional, default true) - if true, standardize data instead of normalization;
##	'Alpha' (optional default 0 for typeT=0, 0 for typeT=1) - the temporal cluster cohesion penalisation factor weight; 
##			-1 only temp component, 0 equal importance, 1 only multidimensional component (Euclidean distance). Proposed in IJAIT 13.
##
## Modified:	
##      13/11/2013 Created this help. Added smooth passage indicator to evaluate typeT=4;
##		03/03/2014 Added the SPassP indicator (penalized with the adjacency matrix);
##		07/03/2014 Returning measures values.
##
## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function [mdvar, tempvar, shannonP, smooth] = evaluateTemporalClustering(varargin)

	[inputfile, clusassign, centroids, ctimestamp, printHeader, outputfid, algoname, algoheader, outputfile, ...
		normalize, stardardize, Alpha, Adj ] = process_options ...
			(varargin, 'inputfile', [], 'clusassign', [], 'centroids', [], 'ctimestamp', [], ...
			'printHeader', false, 'outputfid', [], 'algoname', '', 'algoheader', '', 'outputfile', [],  ...
			'normalize', true, 'stardardize', true, 'Alpha', 0, 'Adj', []);
	
	if (~isempty(outputfile))
		if (~isempty(outputfid))
			warning("Both 'outputfid' and 'outputfile' were specified. Deactivating 'outputfid'");
		endif
		
		# write to filem append to it
		outputfid = fopen(outputfile, 'a');
	endif
	
	if (isempty(outputfid))
		outputfid = 1;
	endif
	
	if (isempty(Adj))
		# initialize Adj to a matrix of one. This will lead to a dummy SPassP = SPass
		Adj = ones(size(centroids, 1));
	endif
	
	if (printHeader)
		fprintf(outputfid, "%s", algoheader);
		fprintf(outputfid, "\tMDvar\tTvar\tSha\tShaP\tSpass\tSPassP\n");
		return;
	endif
						
	if ( isempty(inputfile) || isempty(clusassign) || isempty(centroids) || isempty(ctimestamp) )
		error('Need to give me all the data to work with');
	endif
	
	[files, tags, vals] = loadCsvFileMissing(inputfile, true);
	dtimestamp = vals(:, 1);
	vals(:, 1) = [];
	#if demanded, standardize/normalize data
	if (normalize)
		vals = normalizeMatr(vals, stardardize);
	endif
	noclusters = size(centroids, 1);
	nodata = size(vals, 1);
	datadim = size(vals, 2);
	if (datadim ~= size(centroids, 2))
		error("Data and centroids don't have same number of columns");
	endif
	
	mdvar = 0.0;
	tempvar = 0.0;
	shannon = [];
	shannonP = [];
	
	#calculate the multimensional variance (MDvar)
	datadiff = centroids(clusassign, :);
	datadiff = vals - datadiff;
	datadiff = datadiff .^ 2;
	mdvar = sum(sumskipnan(datadiff)) / nodata;
	
	#calculate the temporal variance (Tvar)
	datadiff = ctimestamp(clusassign);
	datadiff = reshape(dtimestamp, 1, length(dtimestamp)) - datadiff;
	datadiff = datadiff .^ 2;
	tempvar = sum(sumskipnan(datadiff)) / nodata;
	
	#calculate the individual consistency (Shannon entropy) (Sha and ShaP)
	fileU = unique(files);
	for i=1:length(fileU)
		pos = strcmp(files, fileU{i});
		# count in idx how many ovservations are assigned to each individual. Actually, idx contains the position of the last element
		[foo idx] = unique( sort(clusassign(pos)));
		# therefore, we need diferences idx(i)-idx(i-1) to get the count
		probs = diff([0 idx]) / sum(pos);
		shannon = [shannon -sum(log2(probs).*probs)];
		
		nrchanges = length(find(diff(clusassign(pos))));
		
		nrmin = length(probs) - 1;
		shannonP = [shannonP (-sum(log2(probs).*probs))*(1 + (nrchanges-nrmin)/(length(clusassign(pos)) - 1)) ];
	endfor
	shannon = mean(shannon);
	shannonP = mean(shannonP);
	
	#calculate the smooth passage parameter (Spass and SPassP)
	# first, calculate the temporal distance between centroids
	[ctdiff, D] = distanceTemporal('Amultidim', centroids', 'Atimestamp', ctimestamp, 'Bmultidim', centroids', 'Btimestamp', ctimestamp, ...
			'Alpha', Alpha, 'outputfid', 0 );

	smooth = 0;
	smoothP = 0;
	fileU = unique(files);
	for i=1:length(fileU)
		# observations belonging to current individual
		pos = strcmp(files, fileU{i});
		
		# create the vector with the phases through which an entity passes
		# cannot use unique because it orders elements and I lose a->b->a
		idx = find(diff([0 clusassign(pos)]));
		phases = clusassign(pos)(idx);
		
		entsmooth = 0;
		entsmoothP = 0;
		for j=1:length(phases)-1
			entsmooth = entsmooth + ctdiff(phases(j), phases(j+1));
			entsmoothP = entsmoothP + ctdiff(phases(j), phases(j+1)) / Adj(phases(j), phases(j+1));
		endfor
		smooth = smooth + entsmooth / length(phases);
		smoothP = smoothP + entsmoothP / length(phases);
	endfor

	fprintf(outputfid, "%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n", algoname, mdvar, tempvar, shannon, shannonP, smooth, smoothP);
	
	if ( ~isempty(outputfile) )
		fclose(outputfid);
	endif
		
endfunction
