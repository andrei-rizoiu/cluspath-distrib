%% Usage: [lineHeader, columnHeader, values] = loadCsvFile( file, createBin )
%%
%% Function that will read a data file in the CSV text format. The format of the file is
%%
%% <void> 	[tab] colhead1 	[tab] colhead2  [tab] ... [tab] colheadn
%% rowhead1	[tab] val11	[tab] val12	[tab] ... [tab] val1n
%% rowhead2	[tab] val21	[tab] val22	[tab] ... [tab] val2n
%%    ...
%% rowheadm	[tab] valm1	[tab] valm2	[tab] ... [tab] valmn
%%
%% Row and column headers are optional (they can be missing). It is intended to uniformize 
%% the reads from CSV files. If createBin is 1, it creates a binary version of the file 
%% for faster reading a second time. It always tries to read the binary file.
%% The binary file has the name "file.bin" (e.g. tags.csv.bin). file can be a matrix (not a
%% string). In this case it is interpreted as the variable "values" and returned.
%%
%% Modified: 	14/11/2012 Octave 3.6 modified str2double API.
%%		30/11/2012 Octave 3.6 modified modified initialization of cell arrays
%%		25/03/2014 Fixed a bug in determining if the file we receive is a real filename of an octave matrix
%% Author: Marian-Andrei RIZOIU.

function [lineHeader, columnHeader, values] = loadCsvFile( file, createBin, hasLineHeaders )
	
	% need at least one parameter, the file name
	if ( nargin < 1 )
		error("Need to give me a file to read!");
	endif
	
	if ( nargin < 2 )
		createBin = false;	%by default we don't polute the disk with binaries
	endif
	
	if ( nargin < 3 )
		#by default we assume no line headers
		hasLineHeaders = false;	%by default we don't polute the disk with binaries
	endif
	
	% check out if our parameter is a matrix already
	if ( min(min( isprint( file ) )) == 0 || strcmp(typeinfo(file), "string") == 0)	%if at least one character is not readable or the type is not string
		% the parameter is already a matrix, not ment to read it
		columnHeader = [];
		lineHeader = [];
		values = file;
		return;
	endif

	% do we have a binary file we can read?
	binfile = sprintf("%s.bin", file);
	result = dir( binfile );
	if (length(result) == 1 && result.isdir == 0 )
		%only one result and a regular file; lets see what is in it?
		[ vars ] = who("-file", binfile);
		correct = length(vars) == 3;	%need to have the 3 variables here after
		for i=1:length(vars)
			if ( strcmp( vars{i}, "lineHeader" ) == 0 && strcmp( vars{i}, "columnHeader" ) == 0 && strcmp( vars{i}, "values" ) == 0 )
				correct = false;
			endif
		endfor
		
		if ( correct )
			%check that the bin file is not some old version that is not identical to the csv file
			sCSV = stat(file);
			sBIN = stat(binfile);
			
			if (sBIN.mtime > sCSV.mtime)
				%all is good; load the file and exit
				load( binfile );
				return;
			else
				%delete the bin file, since it is obsolete
				unlink( binfile );
			endif
		endif
	endif
	
	%need to check if the file exists
	result = dir(file);
	if (length(result) == 0 || result.isdir == true)
		%before panicking lets see if there is a bzipped version of the file
		bzfile = [ file '.bz2' ];
		bzresult = dir(bzfile);
		if (length(bzresult) ~= 0 && bzresult.isdir == false)
			%bingo: there it is. Let's uncompress it
			bunzip2(bzfile);
		else
			%time to scream
			error( sprintf("File %s cannot be found!", file) );
		endif
	endif
	
	% if here, it means that we need to read the data "a l'ancienne"
	[fid, MSG] = fopen(file, 'r');
	
	if ( fid == -1)
		printf( sprintf("ERROR: Error while loading file! File not found!\n", file) );
		MSG
		columnHeader = [];
		lineHeader = [];
		values = [];
		return;
	endif
	
	lineno = 0;
	line = fgets(fid);
	lineno++;
	%get rid of the end-of-line character
	line( line == "\n") = " ";
	line = strtrim( line );
	columnHeader = strsplit ( line, "\t", true);
	
	% were those that we read really column header or were they already data?
	% for that we test the second element to see if it is a number (as the first might be a line header)
	if ( isnum(columnHeader{2}) > 0 )
		% they were numbers and there is no column header
		columnHeader = [] ;
	else
		% they're alright
		line = fgets(fid);
		lineno++;
	endif
	
	lineHeader = {};
	values = [];
	while ( line ~= -1 )
		testLine = line;
		%read the first field and see if it is a number
		token = testLine(1:index(testLine, "\t")-1);
		if (hasLineHeaders)
			%it is our line header. Add it to the lineHeaders and remove it from input line
			lineHeader = [lineHeader, {token}];
#			testLine = strrep( testLine, [token, "\t"], "")
			testLine = substr(testLine, index(testLine, "\t")+1);
		else 
			if ( !isnum(token) )
				if ( nargin < 3 )
					# if here, it means that we wrongfully assumed that there are no line headersm when the parameter didn't say anything
					# we should restart by assuming line parameters
					[lineHeader, columnHeader, values] = loadCsvFile( file, createBin, true );
					return;
				else
					# if here, the parameters told us that there are no line headers, but there is an error in the data (maybe there are 
					# the headers. Throw an error and die bravely
					error(sprintf("Format error in the given file! The first element on line %d is not a number, but you told me that there are no line headers.", lineno));
				endif
			endif
		endif
		
		%read the rest, which we know now that are numbers
		testLine( testLine == "," ) = ".";		
		val = sscanf(testLine, "%f", [1, inf]);
		
		values = [values; val];
		line =  fgets(fid);
		lineno++;
	endwhile
	fclose(fid);
	
	% if asked we need to save the data we just read in a binary format
	if (createBin == true)
		save("-binary", binfile, "lineHeader", "columnHeader", "values" );
	endif
	
endfunction

function r = isnum(a)
  if ( isnumeric(a) )
    r = 1;
  else
    o = str2double(a);	% TAB  is Ascii 9
    r = !isnan(o);
  endif
endfunction
