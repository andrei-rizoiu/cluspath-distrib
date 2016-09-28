## Usage: extractTemporalComponent(infile, outfile, temptype, args)
##
## Extract the temporal component from the data, using different techniques.
## Parameter:
##	infile - original file
##	outfile	- created file
##	temptype - type of transformation (0 - remove individual average, 1 - differences between observations, 2 - remove individual moving average)
##	args (optional, default [2,2]) - for temptype = 2 how many elements to consider before and after the current one in the average.

function extractTemporalComponent(infile, outfile, temptype, args)
	
	if (nargin < 3)
		temptype = 0;
	endif
	
	[f,t,v] = loadCsvFileMissing(infile, true);
	newv = zeros(size(v));
	
	if (temptype == 0)
	    # separate the year from the data
		temp = v(:, 1);
		v(:, 1) = [];
		newv = zeros(size(v));
		
		# list of entities
		countr = unique(f);
		for i=1:length(countr)
			pos = strcmp(f, countr{i});

			[ssum, ni] = sumskipnan(v(pos, :));
			ni( ni == 0 ) = eps;
			avr = ssum ./ ni;
			
			newv(pos, :) += repmat(avr, sum(pos), 1);
		endfor
		
		newv = v - newv;
		newv = [ temp newv];
	endif
	
	if (temptype == 1)
		temp = v(:, 1);
		v(:, 1) = [];
		newv = zeros(size(v));
		
		father = zeros(1, columns(v));	#the non-NaN element position before the current one
		for i=1:rows(v)
			for j=1:columns(v)
				#current position calculation
				if ( father(j) == 0 || strcmp(f{i}, f{father(j)}) == 0 )
					newv(i,j) = NaN;
				else
#					printf('v(%d, %d)=%f; ', i,j,v(i,j));
#					printf('father(%d)=%d; ', j,father(j));
#					printf('temp(%d)=%d; ', i, temp(i));
#					printf('temp(father(%d))=%d\n', j, temp(father(j)));					
					newv(i, j) = (v(i, j) - v(father(j), j) ) / (temp(i) - temp(father(j)) );
				endif
				
				#father update
				if ( ~isnan(v(i,j) ) )
					#if current element is non NaN, then it will be the father of next ones
					father(j) = i;
				endif
			endfor
		endfor

		newv = [ temp newv];
	endif
	
	if (temptype == 2)
	
		temp = v(:, 1);
		v(:, 1) = [];
		newv = zeros(size(v));
	
		if (nargin < 4)	#mobile average of window 5 by default
			before = 2;
			after = 2;
		else
			before = args(1);
			after = args(2);
		endif
		
		for i=1:rows(v)
			index = [i-before:i+after];
			
			#remove aberant indexes
			index(index < 1) = [];	
			index(index > rows(v)) = [];
			
			#verify that it belongs to the same country
			mycountr = f{i};
			index(strcmp(mycountr, f(index)) == 0) = [];
			
			#calculate the moving average
			[ssum, ni] = sumskipnan(v(index, :));
			ni( ni == 0 ) = eps;
			avr = ssum ./ ni;
			
			newv(i, :) = avr;
		endfor
		
		newv = v - newv;
		newv = [ temp newv];
	endif
	
	saveCsvFile(f, t, newv, outfile);
endfunction
