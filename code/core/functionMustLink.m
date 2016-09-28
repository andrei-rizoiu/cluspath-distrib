## Usage: [MLP] = functionMustLink(tdata, files, clusassign, noclusters, args)
##
## Function used to calculate the contiguity penalisation component in the objective function. Proposed in ICTAI 2012.
## It calculats the penality based on semi-supervised mustlink constraints: mustlink constraints are set between all observations belonging to the same entity.
## If two of these observations are partitioned in different clusters, a penalty is inflicted, inversly proportional to the time difference betweent the two.
## 
## The function is (chosen by 'normal_distrib_function' parameter): 
##		Beta * min(1, 1 / |t(di) - t(dj)| )
## or 
##		Beta * exp( -0.5 .* ((x / Delta) .^ 2) ) - normal distribution function.
## Complete function: 
##		Beta * sum(dk) [min(1, 1 / |t(di) - t(dj)| ) * (clus(di) == clus(dk)) ]
##
## Outputs: MLP(i,j) is the penality associated with the element di clustered in cluster cj
##
## Parameters: 
##		'tdata' (default []) - timestamp of observations (dtimestamp);
##		'files' (default []) - files{i} is the entity to which observation di is associated;
##		'clusassign' (default []) - the cluster assignement of observations at previous iteration;
##		'noclusters' (default []) - number of clusters;
##		'Beta' (default 1) - the temporal entity cohesion penalisation factor weight. Set to zero for no cohesion information.
##		'Delta' (default 3) - the width of the contiguous segmentation function. Set Delta = 3 in order to have around 0.4Beta at a DeltaT=4;
##		'typeT' (detault 1) - see defition in "temporalKMeans.m"
##		'normal_distrib_function' (default true) - true to use the normal distribution-inspired function above
##					set false to use Beta * min(1, 1 / |t(di) - t(dj)| )
#3		'Adj' (default []) - the adjacency matrix. Needed only for typeT=5.
##
## Modified:	29/09/2013 Added this help, adapted for typeT = 3 (graph structure inside clustering);
##		01/11/2013 Made our proposal by default, should not need more intervention for future typeT
##		27/11/2013 Added support for typeT = 5
##
## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function [MLP] = functionMustLink(varargin)

	[tdata, files, clusassign, noclusters, Beta, Delta, typeT, normal_distrib_function, Adj ] = process_options (varargin, ...
		'tdata', [], 'files', [], 'clusassign', [], 'noclusters', [], 'Beta', 1, 'Delta', 3, ...
		'typeT', 1, 'normal_distrib_function', true, 'Adj', []);

	
	if ( isempty(tdata) || isempty(files) || isempty(clusassign) || isempty(noclusters) )
		error('You did not give me enought data to work with');
	endif
	
	if ( isempty(Adj) && (typeT == 5) ) 
		error(sprintf("The selected algorithm (typeT=%d) requires setting the adjacency matrix!", typeT));
	endif
	
	if ( (typeT == 5) && ~normal_distrib_function )
		warning(sprintf("The selected algorithm (typeT=%d) requires the normal distribution-inspired penalty function.Setting this parameter!", typeT));
		normal_distrib_function = true;
	endif
	
	######################### END consistency check ###############################
	
	noobservations = length(tdata);
	MLP = zeros(noobservations, noclusters);
	if (Beta == 0)
		return;
	endif
	
	## to calculate MLP(i,j), we must find the entities dk so that clus(di) != clus(dk) and ind(di) == ind(dk)
	
	for i=1:noobservations
		slct = strcmp(files, files{i});	#which observations correspond to the same entity as i ind(di) == ind(dk)
		slct = repmat(slct', 1, noclusters);

		## at this end of the calculation, temp(p,q) represents the temporal penalty inflicted by the pair of observations
		## (x_i, x_p), if the observation x_i is assigned to the cluster p. The assignement for the observation x_q
		## remains as in the previous iteration. If x_i and x_p do not share the same entity (info in slct), then
		## temp(p, :) = NaN. OBS: a pair is taken into account twice, since the direction is not considered. We calculate
		## the penalty from an observation to any observation sharing the same entity (either before or after the current 
		## observation).
		temp = repmat(vec(clusassign), 1, noclusters);
		temp(i, :) = [1:noclusters];				#possible assignement of observation i
		temp = (temp ~= repmat(temp(i, :), noobservations, 1)) + 1 - 1;	#wheather elements are in the same or different cluster
		temp = temp .* slct;	#wheather the observations correspond to the same entity
		temp( temp == 0) = NaN;
		temp = temp .* repmat(vec(tdata), 1, noclusters);
		#pe caiet, temp este X (temp(j, k) is the vector time difference between obs i and j when i is put in cluster k - if i and j in difference clusters and same entity)
		temp = tdata(i) - temp;	
		
		switch typeT
			case {2}
				# the windows based function proposed in LIN06
				temp = abs(temp);
				temp = (temp <= Delta) + 1 - 1;	#just set it to 1 is the timedifference less than a window
				
			case {5}
				# we are setting forward precedence links, from the current observation (x_i) to all the other observations (x_p)
				# We will, afterwards, filter the following observations:
				#	- which do not belong to the same entity as x_i (set to NaN in temp using the "slct" variable, so no need here);
				#	- which are in the same cluster as x_i (also set to NaN in temp, above);
				#	- which are antecedents of x_i (will filter them later).
				#
				# The forward links is done between the new cluster of x_i and the cluster of x_p, therefore on the link in the adjacency matrix
				# given by "q" and clusassign(p), where q takes all the values from 1 to noclusters (columns in temp).
				# we need to weight the temp(p, q) with the term 1 - Adj(q, clusassign(p))
				# q is the new assignement of the observation x_i and clusassign(p) is the old
				# cluster assignement of observation x_p, if tdata(i) < tdata(p).
				adj_pen = ones(noobservations, noclusters);

				# Adj(q, clusassign(p)) = Adj(q, clusassign) {since P = 1:noobservations} = Adj(:, clusassign) {since q = 1:noclusters}
				adj_pen = 1 - Adj(:, clusassign)' .^ 2;

				# only interested in forward must-links from x_i to x_p, where x_p^t>x_i^t (x_p which are succesors of x_i)
				# therefore only those "temp" < 0
				adj_pen(temp >= 0) = NaN;
				
				# no more filtering is required, since verification for being in the same entity and different clusters
				# has been already done in the temp. Elements for which this condition is not fulfilled are NaN
				# first calculate the rest of the penalty factor, just in other typeT's (in otherwise branch)
				# we force the normal distribution-inspired penalty function
				temp = exp( -0.5 .* ((temp / Delta) .^ 2) );
				# simply apply the penalty				
				temp = temp .* adj_pen;

			otherwise
				# default function, proposed by us in ICTAI12
				temp = abs(temp);
				if (normal_distrib_function)
					temp = exp( -0.5 .* ((temp / Delta) .^ 2) );
				else
					temp( temp < 1) = 1;
					temp = 1 ./ temp;
				endif
		endswitch
		
		
		[temp, Ni] = sumskipnan(temp);
		MLP(i, :) = temp;
	endfor

	MLP = MLP * Beta;
	
endfunction
