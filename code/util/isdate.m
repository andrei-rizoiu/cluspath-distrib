## Usage: dater = isdate(a)
##
## Function takes a string or cell array of strings as parameter 
## and returns an array of "datenum" (see datenum) date 
## representation if it is a date or NaN if it is not a date.
## Tries iteratively several date formats.
##
## Created: 12/11/2012

## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function dater = isdate(a)
	if ( iscell(a) )
		dater = [];
		for i = 1:length(a)
			res = isdate_single(a{i});
			dater = [ dater res ];
		endfor
	else
		dater = isdate_single(a);
	endif	
endfunction

# function to do it for only one string.
function dater = isdate_single(a)
	# put only current format for speedup
#	formats = {"yyyy-mm-dd HH:MM:SS"};
	formats = {"yyyy-mm-dd HH:MM:SS" "yyyy-mm-dd"};
	for i = 1:length(formats)
		try
			dater = datenum(datevec(a, formats{i}));
		catch
			dater = NaN;
		end_try_catch
		
		#found a format that matches, just exit
		if ( !isnan(dater) )
			return
		endif
	endfor
endfunction
