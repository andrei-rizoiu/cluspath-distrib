## Usage: [f, t, v] = addClusterInfo(infile, info, outfile)
##
## Adds the cluster info to the file, so that it can be visualized in Weka.

function [f, t, v] = addClusterInfo(infile, cinfo, outfile)

	[f, t, v] = loadCsvFileMissing(infile);
	v = [v cinfo'];
	t = [t {'ClusterOctave'}];
	
	saveCsvFile (f, t, v, outfile);
endfunction
