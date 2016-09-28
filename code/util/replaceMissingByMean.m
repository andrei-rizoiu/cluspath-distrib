## Usage: out = replaceMissingByMean(in)
##
## Replaces missing valeus (NaN in octave) by the mean(average) of the other values present.

## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function data = replaceMissingByMean(data)

	[averData, yData] = sumskipnan(data);
	averData = averData ./ yData;
	
	repl = isnan(data);
	multip = repmat(averData, rows(data), 1);
	repl = repl .* multip;
	
	data(isnan(data) == 1) = 0;
	data = data + repl;

endfunction
