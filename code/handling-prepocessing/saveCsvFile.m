## Usage: saveCsvFile( lineHeader, columnHeader, vals, file, columnHeadersAtTheEnd, sep )
##
## Function that will save in a data file in the text format.
## Note: line and column headers can be absent
##
## Parameters:
##	lineHeader (optional) - the line headers;
##	columnHeader (optional) - the column headers;
##	vals - the matrix of values to write in the file;
##	file - the name of the file to write to;
##	columnHeadersAtTheEnd (default false) - put column headers at the end of line, instead of beginning.
##	sep (default '\t') - field separation character
##
## Modified:	04/11/2013 - added this help, parameter check, column headers at the end.
##
## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function saveCsvFile( lineHeader, columnHeader, vals, file, columnHeadersAtTheEnd, sep )

	if (isempty(vals) || isempty(file))
		error('The vals and file parameters are compulsory!');
	endif
	
	if (nargin < 5) 
		columnHeadersAtTheEnd = false;
	endif
	
	if (nargin < 6) 
		sep = "\t";
	endif

	[fid, MSG] = fopen(file, 'w');
	
	if ( fid == -1)
		MSG
		return;
	endif
	
	if (columnHeadersAtTheEnd && ( columns(columnHeader) > columns(vals) ))
		% we need to take the first column header and put it at the end
		header = columnHeader{1};
		columnHeader(1) = [];
		columnHeader = [columnHeader {header} ];
	endif
	
	% write the column headers, if they exist
	if (columns(columnHeader) > 0 )
		if ( columns(columnHeader) == columns(vals) )
			fprintf(fid, "%s", sep);
		endif
		for i=1:columns(columnHeader)-1
			fprintf(fid, "%s%s", columnHeader{i}, sep);
		endfor
		fprintf(fid, "%s\n", columnHeader{columns(columnHeader)});
	endif
	
	%write each line
	for i=1:rows(vals)
		if ( ( columns(lineHeader) > 0 ) && !columnHeadersAtTheEnd)
			fprintf(fid, "%s%s", lineHeader{i}, sep);
		endif
		for j=1:columns(vals)-1
			fprintf(fid, "%.10f%s", vals(i, j), sep);
		endfor
		fprintf(fid, "%.10f", vals(i, end ));
		
		if ( ( columns(lineHeader) > 0 ) && columnHeadersAtTheEnd)
			fprintf(fid, "%s%s", sep, lineHeader{i});
		endif
		
		fprintf(fid, "\n");
	endfor
	
	fclose(fid);
end
