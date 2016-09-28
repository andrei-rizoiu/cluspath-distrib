## Usage: cumulativeBinary(filename, outfile)
##
## Given a file in the binary format (which was created with Weka
## and reput into form using "wekacsv2octave.sh" script), this turns
## it into the cumulative bunary format. Ex:
## a=0, a=1 and a=2 becomes
## a<=0, a<=1 and a<=2
##
## Matrix:
## 1 0 0	1 1 1
## 1 0 0	1 1 1
## 0 1 0	0 1 1
## 0 1 0	0 1 1
## 0 0 1	0 0 1
## 0 0 1	0 0 1

## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function cumulativeBinary(filename, outfile)
	if (nargin < 2)
		usage('cumulativeBinary(filename, outfile)');
	endif
	
	[f,t,v] = loadCsvFile(filename);
	
	oldTag = [];
	tag = [];
	for i=1:length(t)
		tag = t{i};
		
		#get the initial tag name
		[tagIni, tagRem] = strtok(tag, "=");
		[oldTagIni, oldTagRem] = strtok(oldTag, "=");
		
		if ( strcmp(tagIni, oldTagIni) == 1 )
			#we are dealing with the same initial tag with two consecutive values
			[newVal, foo] = strtok(tagRem, "=");
			newName = sprintf("%s<=%s", tagIni, newVal);
			t{i} = newName;
			
			v(:, i) = v(:, i) + v(:, i-1);
		endif
		
		oldTag = tag;
	endfor
	
	saveFile(f, t, v, outfile);
endfunction
