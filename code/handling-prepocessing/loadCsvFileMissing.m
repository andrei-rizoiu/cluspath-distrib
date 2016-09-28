## Usage: [lineHeader, columnHeader, values] = loadCsvFile( file, createBin )
##
## Function that will read a data file in the CSV text format. The format of the file is
##
## <void> 	[tab] colhead1 	[tab] colhead2  [tab] ... [tab] colheadn
## rowhead1	[tab] val11	[tab] val12	[tab] ... [tab] val1n
## rowhead2	[tab] val21	[tab] val22	[tab] ... [tab] val2n
##    ...
## rowheadm	[tab] valm1	[tab] valm2	[tab] ... [tab] valmn
##
## Row and column headers are optional (they can be missing). It is intended to uniformize 
## the reads from CSV files. If createBin is 1, it creates a binary version of the file 
## for faster reading a second time. It always tries to read the binary file.
## The binary file has the name "file.bin" (e.g. tags.csv.bin). file can be a matrix (not a
## string). In this case it is interpreted as the variable "values" and returned.
##
## This function takes care of missing values. Puts "?" (missing values) as NaN.
## Warning: use this only if needed, since we lost vector operations and there will be 
## performance penalties.
##
## Also tries to interpret values that cannot be transformed into numbers as data.
## It uses the "isdate" function that tries iteratively different date formats. This
## approach is useful when the date is in a date format (as the Technicolor dataset).
##
## OBSERVATION: the file parameter can also be directly the matrix to be loaded. In this case,
##	we return [ [], [], file]. This is a bogus behavious in order to allow seamless reading
##	from files and from structures already in memory.
##
## Modified: 	12/11/2012 added date format support; Octave 3.6 modified str2double API.
##		12/11/2013 added support for Weka generated files, which have "," as field separator
##		04/03/2014 corrected bug when the file was a matrix and contained NaNs
##
## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function [lineHeader, columnHeader, values] = loadCsvFileMissing( file, createBin )

	DEBUG = 0;
	
	# need at least one parameter, the file name
	if ( nargin < 1 )
		error("Need to give me a file to read!");
	endif
	
	if ( nargin < 2 )
		createBin = false;	#by default we don't polute the disk with binaries
	endif
	
	# check out if our parameter is a matrix already
	if ( sum(sum(isnan(file))) > 0 || rows(file) > 1 || min(min( isprint( file ) )) == 0 )	#if at least one character is not readable
		% the parameter is already a matrix, not ment to read it
		columnHeader = [];
		lineHeader = [];
		values = file;
		return;
	endif

	# do we have a binary file that we can read?
	binfile = sprintf("%s.bin", file);
	result = dir( binfile );
	if (length(result) == 1 && result.isdir == 0 )
		#only one result and a regular file; lets see what is in it?
		[ vars ] = who("-file", binfile);
		correct = length(vars) == 3;	#need to have the 3 variables here after
		for i=1:length(vars)
			if ( strcmp( vars{i}, "lineHeader" ) == 0 && strcmp( vars{i}, "columnHeader" ) == 0 && strcmp( vars{i}, "values" ) == 0 )
				correct = false;
			endif
		endfor
		
		if ( correct )
			#check that the bin file is not some old version that is not identical to the csv file
			sCSV = stat(file);
			sBIN = stat(binfile);
			
			if (sBIN.mtime > sCSV.mtime)
				#all is good; load the file and exit
				load( binfile );
				return;
			else
				#delete the bin file, since it is obsolete
				unlink( binfile );
			endif
		endif
	endif
	
	#need to check if the file exists
	result = dir(file);
	if (length(result) == 0 || result.isdir == true)
		#before panicking lets see if there is a bzipped version of the file
		bzfile = [ file '.bz2' ];
		bzresult = dir(bzfile);
		if (length(bzresult) ~= 0 && bzresult.isdir == false)
			#bingo: there it is. Let's uncompress it
			bunzip2(bzfile);
		else
			#time to scream
			error( sprintf("File %s cannot be found!", file) );
		endif
	endif
	
	# if here, it means that we need to read the data "a l'ancienne"
	[fid, MSG] = fopen(file, 'r');
	
	if ( fid == -1)
		printf( sprintf("ERROR: Error while loading file! File not found!\n", file) );
		MSG
		columnHeader = [];
		lineHeader = [];
		values = [];
		return;
	endif
	
	line = fgets(fid);
	%get rid of the end-of-line character
	line( line == "\n") = " ";
	line = strtrim( line );
	columnHeader = strsplit ( line, "\t", true);
	replaceCommaByTab = 0;
	if (length(columnHeader) == 1)
		# for Weka originating CSVs, separator is ",". Replace it
		line( line == ",") = "\t";
		replaceCommaByTab = 1;
		fprintf(DEBUG, "Field separator detected as being \",\" instead of tabulator. Replacing...\n");
		columnHeader = strsplit ( line, "\t", true);
	endif
	
	% were those that we read really column header or were they already data?
	% for that we test the second element to see if it is a number (as the first might be a line header)
	if ( isnum(columnHeader{2}) > 0 )
		% they were numbers and there is no column header
		columnHeader = [] ;
	else
		% they're alright
		line = fgets(fid);
	endif
	
	lineHeader = {};
	values = [];
	linecount = 1;
	fprintf(DEBUG, "Processing line:\t");
	while ( line ~= -1 )

		linecount++;
		if ( mod(linecount, 100) == 0 )
			fprintf(DEBUG, "%d\t", linecount);
		endif
		
		if (replaceCommaByTab)
			line( line == ",") = "\t";
		endif
		testLine = line;
		#read the first field and see if it is a number
		file = testLine(1:index(testLine, "\t")-1);
		if ( !isnum(file) )
			#nope, it is our line header. Add it to the lineHeaders and remove it from input line
			lineHeader = [lineHeader, {file}];
			testLine = strrep( testLine, [file, "\t"], "");
		endif
		
		# read the rest, which we know now that are numbers
		# the next line is to take into account decimal separator for French and Romanian languages, which is comma
		testLine( testLine == "," ) = ".";	
		numbers = strsplit ( testLine, "\t", true);	

		# we start the magic: convert to double the strings
		val = str2double(numbers);
		# where convertion issued a NaN, try to treat it as date and convert the date to strings
		datepos = isnan(val);
		dates = isdate(numbers( datepos ));
		[res, nval] = sumskipnan( [val(datepos) ; dates] );
		res (nval == 0) = NaN;
		val(datepos) = res;
#		val = nansum( [str2double(numbers) ; isdate(numbers)] );
		
		values = [values; val];
		line =  fgets(fid);
	endwhile
	fclose(fid);
	
	# if asked we need to save the data we just read in a binary format
	if (createBin == true)
		save("-binary", binfile, "lineHeader", "columnHeader", "values" );
	endif
	
	fprintf(DEBUG, "\n");
	
endfunction

# checks if a is a valid number
function r = isnum(a)
  if ( isnumeric(a) )
    r = 1;
  else
    o = str2double(a);
    r = !isnan(o);
  endif
endfunction

