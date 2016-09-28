## Usage: [interphi, counts, ctimestamp, arcs, entities, startsat] = interPhi('files', files, 'dtimestamp', dtimestamp, 'clusassign', clusassign, 'filtcycle', filtcycle, 'ctimestamp', ctimestamp);
##
## Based on the entity name, dtimestamp and cluster assignement for each observation, this function
## calculates the inter_PHI function between clusters/phases. The inter_PHI(\mu_p, \mu_q) returns
## the number of entities which have observations in both clusters \mu_p and \mu_q and those in
## cluster \mu_p precede temporally those in \mu_q
##
## Parameters: 	'filename' or ('files', 'dtimestamp', 'clusassign') the data needed
##		'ctimestamp' (default []) - the timestamps of centroids. If not timestamp is given (default []), then recompute
##					timestamps by averaging the timestamps of the observations in clusters.
##		'filtcycle' (default false) controls wheather cycles 1 -> 2 -> 1 are moved
##
## It outputs several values:
##	- interphi - the interphi matrix, where interphi(i, j) is the intersection distance between clusters i and j.
##	- counts - the counts matrix, counts(i,j) is the number of entities passing from cluster i to cluster j
##	- ctimestamp - recalculates the timestamps of centroids, only if they were not given by parameters;
##	- arcs - arcs{i, j} is a list of entities that pass from state i to state j
##	- entities - are the names of entities, in ascending order;
##	- starts - a vector of the same dimension as the list of entities, so that startsat(i) show the first
##		cluster for entity i.
##
## Modified:	16/11/2012 Revisited, added this help
## 	   	22/04/2013 Revisited, added the possibility to choose if to filter cycles; changed name
##		29/09/2013 Added support for receiving the ctimestamp as parameter and calculating it only if necesary.
##				Parameter passing using "process_options" function, improved help.

## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function [interphi, counts, ctimestamp, arcs, entities, startsat] = interPhi(varargin)

	[filename, files, dtimestamp, clusassign, ctimestamp, filtcycle] = ...
		process_options(varargin, 'filename', [], 'files', [], 'dtimestamp', [], 'clusassign', [], ...
		'ctimestamp', [], 'filtcycle', false );
		
	if ( ~isempty(filename) )
		[files, tags, vals] = loadCsvFileMissing(filename, true);
		dtimestamp = vals(:, 1);
		clusassign = vals(:, end);
		vals = [];
	endif
	
	if ( isempty(files) || isempty(dtimestamp) || isempty(clusassign) )
		error("I don't have all the data I need to work!");
	endif
	
	clust = unique(clusassign);
	counts = zeros( length(clust), length(clust) );
	interphi = zeros( length(clust), length(clust) );
	arcs = cell(size(counts));
	
	######################################## ENDED INITIALIZATION ######################
	
	# if the ctimestamp was not given, calculate ctimestamp by averaging the timestamps of the observations in clusters.
	if ( isempty(ctimestamp) )
		ctimestamp = zeros(1, length(clust));
		for i = 1:length(clust)
			pos = find(clusassign == i);
			ctimestamp(i) = sum(dtimestamp(pos)) / length(pos);
		endfor
	endif
	
	entities = unique(files);
	startsat = zeros(1, length(entities));
	for i = 1:length(entities)
		pos = find(strcmp(entities{i}, files));	# shows which observations belong to current entity
		centr = clusassign(pos)';
		node1 = centr(1);
		
		# eliminate doubles - OBS: cannot use unique as it return elements in ascending order
		# I need to obtain from 1112221113333311122244422 the chain 12131242
		centr_new = centr(1);
		for j = 2:length(centr)
			if (centr(j) ~= centr(j-1))
				centr_new = [centr_new, centr(j)];
			endif
		endfor
		centr = centr_new;
		
		if (filtcycle)
			#eliminate cycles "a" -> "b" -> "a", by removing "b" and obtaining only "a"
			centr_new = [];
			adding = [];
			ads = "";
			cod = 0;
			for j = 1:length(centr)
				cod = cod - 1;
				if (cod < 0)
					if (j <= length(centr)-2) 
						if (centr(j) ~= centr(j+2))
							centr_new = [centr_new, centr(j)];
							adding = [adding, {ads}];
							ads = "";
						else
							cod = 1; #skip the next one
							ads = sprintf("(a/r %d)", centr(j+1));
						endif
					else
						#there can be no cycle now
						centr_new = [centr_new, centr(j)];
						adding = [adding, {ads}];
					endif
				endif
			endfor
			centr = centr_new;
		endif
		
		startsat(i) = node1;

		for j = 2:length(centr)
			node2 = centr(j);
			interphi(node1, node2) = interphi(node1, node2) + 1;			

			# i need to verify that arcs{node1, node2} does not already contain my entity, given that cycles
			# are not necessarily omitted. It is due to evolutions "a -> b -> a -> c", so the entity appears twice
			# for "a"
			if ( sum(centr(2:j-1) == node2 ) == 0 ) 	
				counts(node1, node2) = counts(node1, node2) + 1;
				extra = "";
				if (length(arcs{node1, node2}) ~= 0)
					extra = ", ";
				endif
				arcs{node1, node2} = [arcs{node1, node2}, extra, entities{i} ];
			endif

			node1 = node2;
		endfor
	endfor

	# calculate the interphi function which is 1 - counts/number-of-entities
	normal = size(entities)(2) ; # normalization factor
	interphi = 1 - (counts / normal); # OBS: the counts variable earlier is the official counts, where a country does not
						# appear twice, regardless of the value for filtCycle. This is because counts is printed
						# in the "generateGraphFile" function
	
endfunction
