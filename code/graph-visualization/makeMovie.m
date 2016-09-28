function makeMovie(Adj, ctimestamp, inputfile)

	folder = [ tempdir() "/octave-movie/" ] ;
	mkdir(folder);
	no_arcs = length(find(Adj)) ;
	no_arcs = min(no_arcs, 100) ;
	
	for i = 0:no_arcs
		filename = sprintf("%s/output-graph-%05darcs.dot", folder, i);
		generateGraphFile('inputfile', inputfile, 'outputfile', filename, 'edgedescription', false, 'Adj', Adj, 'ctimestamp', ctimestamp, 'no_arcs', i, 'removeIsolatedNodes', false);
	endfor
endfunction
