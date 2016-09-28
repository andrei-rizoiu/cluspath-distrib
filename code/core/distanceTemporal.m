## Usage: [dTemp, d] = distanceTemporal('Amultidim', [], 'Atimestamp', [], 'Bmultidim', [], 'Btimestamp', [], 'files', [], 'clusassign', [], 'Alpha', [], 'Beta', [], 'Delta', [], 'missing', false, 'typeT', 1, 'outputfid', 1, 'deltaXmax', -1, 'deltaTmax', -1 )
##
## Calculates the Temporal Aware distance, that takes into account the temporal and the multidimensional components. Published in ICTAI12.
## The temporal aware distance calculated here can contain the congiguity penalties if "clusassign" is not empty. See publication for details.
## Returns both the temporal aware distance including the temporal penalities and the simple euclidian distance:
##	-> dTemp - dTemp(i,j) is the temporal aware distance between individual i in the set A and individual j from set B (including contiguity penalties)
##	-> d - d(i,j) is the SQUARE of the Euclidean distance between individual i in the set A and individual j from set B
##
## Parameters:
##	- Amultidim - matrix containing individuals in columns
##	- Atimestamp - unidimensional vectors of timestaps for the Amultidim
## 	- Bmultidim - matrix containing individuals in set B in columns
##	- Btimestamp - unidimensional vectors of timestaps for the individuals in set B
##	- files - the individuals to which observations belong (needed only for contiguity component)
##	- clusassign - cluster assignement at the previous iteration (clusassign(i) is the cluster no of element i) (needed only for contiguity component)
##	- Alpha (optional default 0) - the temporal cluster cohesion penalisation factor weight. 
##			-1 only temp component, 0 equal importance, 1 only multidimensional component. 
##	- Beta (optional default 0.003) - the temporal individual cohesion penalisation factor weight. Zero for no temporal information.
##	- Delta (optional, default 3) - the width of the contiguous segmentation function;
##	- missing (optional, default false) - does your data contain missing values? If missing data is detected, that it is set automatically to true.
##	- typeT (default 1) - type of temporal function used (see definition in "temporalKMeans.m" file)
##	- outputfid (optional, default 1) - the FID to which to print debug messages. 1 is the normal output, 2 is the error, 0 is silent (no debug)
##	- deltaXmax (optional, default -1) - the spatial diameter. -1 means no diameter, the max in the matrix is taken
##	- deltaTmax (optional, default -1) - the temporal diameter. -1 means no diameter, the max in the matrix is taken
##	- Adj (optional, default []) - the adjacency matrix, needed for typeT=5. Passed to the "functionMustLink" function
##
## Examples: 
## 	[dTemp, d] = distanceTemporal('Amultidim', Amultidim', 'Atimestamp', Atimestamp, 'Bmultidim', Bmultidim', 'Btimestamp', Btimestamp); - calculates temporal-aware 
##				distance between A and B (with their respective 2 components). All parameters are to default values.
##
## Modified:	27/09/2013 Improved the help.
##		29/09/2013 Added support for typeT=3.
##		01/11/2013 Added support for generic temporal-aware distance (initially only between individuals and centroids). Needed for distance between centroids.
##		27/11/2013 Added support for typeT=5.
##
## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>


function [dTemp, d] = distanceTemporal(varargin)

	[Amultidim, Atimestamp, Bmultidim, Btimestamp, files, clusassign, Alpha, Beta, Delta, missing, typeT, outputfid, deltaXmax, ... 
	deltaTmax, Adj] = ...
		process_options (varargin, 'Amultidim', [], 'Atimestamp', [], 'Bmultidim', [], 'Btimestamp', [], 'files', [], ... 
		'clusassign', [], 'Alpha', 0, 'Beta', 0.003, 'Delta', 3, 'missing', false, 'typeT', 4, 'outputfid', 1, ...
		'deltaXmax', -1, 'deltaTmax', -1, 'Adj', [] );

	################ Sanity checks ##############################

	if ( isempty(Amultidim) || isempty(Atimestamp) || isempty(Bmultidim) || isempty(Btimestamp) )
		error('You did not give me enought data to work with');
	endif
	
	if (~isempty(clusassign) && isempty(files) )
		error('You need to give me the file information if I am to calculate the congiguity component');
	endif
	
	# Atimestamp and Btimestamp need to be unidimensional
	if (sum(size(Atimestamp) == 1) < 1)
		error('Atimestamp needs to be unidimensional matrix');
	endif
	
	if (sum(size(Btimestamp) == 1) < 1)
		error('Btimestamp needs to be unidimensional matrix');
	endif
	
	# we want Atimestamp as a column vector, whereas Btimestamp as a row vector
	# this is due to legacy programming from when the calculation was between individuals and centroids
	Atimestamp = vec(Atimestamp);
	Btimestamp = vec(Btimestamp)';
	
	# Amultidim and Atimestamp need to have the same number of individuals (individuals in columns)
	if ( size(Amultidim, 2) ~= length(Atimestamp) )
		error('Amultidim and Atimestamp need to have the same number of individuals (remember, individuals in columns in Atimestamp)');
	endif
	
	# Bmultidim and Btimestamp need to have the same number of individuals (individuals in columns)
	if ( size(Bmultidim, 2) ~= length(Btimestamp) )
		error('Bmultidim and Btimestamp need to have the same number of individuals (remember, individuals in columns in Btimestamp)');
	endif
	
	# finally, Amultidim and Bmultidim need to be described in the same space (same number of features)
	if ( size(Amultidim, 1) ~= size(Bmultidim, 1) )
		error('Amultidim and Bmultidim need to have the same number of features (remember, individuals in columns in Btimestamp)');
	endif
	
	################ END of sanity checks #######################

	# calculate gammaD and gammaT
	gammaD = 1 + Alpha ; 
	if (gammaD > 1)
		gammaD = 1;
	endif

	gammaT = 1 - Alpha ; 
	if (gammaT > 1)
		gammaT = 1;
	endif

	# first calculate the Euclidean distances between elements in of dataset and individuals in set B
	if ( (sum([1 3 4] == typeT) > 0) && (gammaD == 0))
		# in the case of typeT 1, 3, 4, there is no need for any calculation
		# since it will be destroyed later because gammaD is zero, simply create a zero matrix
		d = zeros(columns(Amultidim), columns(Bmultidim));
	else
		if ( ! missing )
			d = distance(Amultidim, Bmultidim);
			if ( sum(sum(isnan(d))) > 0)
				missing = true;
				#there are NaNs in the distance, maybe someone forgot to tell us that we have missing data
				d = distanceMissing(Amultidim, Bmultidim);
			endif
		else
			d = distanceMissing(Amultidim, Bmultidim);
		endif
		d = d .^ 2;
	endif
	
	## d(i,j) is the square euclidian distance between Amultidim(:, i) and centroid(:, j)
	
	# do extra calculations on the Euclidean distance, if needed
	switch typeT
		case {0, 2}
			dTemp = d;
	
		case {1, 3, 4, 5}
			fcts = ones(size(d));
			if (gammaT ~= 0)
				x1 = repmat(Atimestamp, 1, columns(Bmultidim));
				x2 = repmat([1:columns(Bmultidim)], rows(Atimestamp), 1);
				dts = x1 - Btimestamp(x2);
				fcts = dts .^ 2;
				if (deltaTmax ~= -1)
					mx = deltaTmax ^ 2;
				else
					mx = max(max(fcts));
				endif
				fcts = 1 - gammaT * (fcts / (mx));
			endif

			if (deltaXmax ~= -1)
				mx = deltaXmax ^ 2;
			else
				mx = max(max(d));
			endif
			#what if gammaD is zero and d is only zeros?
			if (mx == 0)
				mx = eps;
			endif
			dm = 1 - gammaD * (d / mx);

			dTemp = 1 - dm .* fcts;
		
		otherwise
			error(sprintf("I don't recognize 'typeT=%d' Did you give me a wrong typeT algorithm to work with???", typeT));
	endswitch
	
	# calculate the "mustlink" component
	MLP = zeros(size(d));
	if ( ~isempty(clusassign) ) #it is always empty for the first iteration, when no membership exists already
		[MLP] = functionMustLink('tdata', Atimestamp, 'files', files, 'clusassign', clusassign, ...
				'noclusters', columns(Bmultidim), 'Beta', Beta, 'Delta', Delta, 'typeT', typeT, 'Adj', Adj );
	endif
	
	# DEPRECATED UNPUBLISHED
	# calculate the "cannotlink" component, the WITHIN component which handles cluster cohesion
	# only for typeT = 0
	CLP = zeros(size(d));
	if ( ~isempty(clusassign) && (typeT == 0) ) #it is always empty for the first iteration, when so membership exists already
		[CLP] = functionCannotLink(Atimestamp, clusassign, columns(Bmultidim), Alpha);
	endif

	# add the temporal components
	fprintf(outputfid, "Means: TA distance: %.3f, contiguity: %.3f\n", mean(mean(dTemp)), mean(mean(MLP)));
	dTemp = dTemp + MLP + CLP;

endfunction

