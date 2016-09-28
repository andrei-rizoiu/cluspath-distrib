## Author: Marian-Andrei RIZOIU <andrei.rizoiu@gmail.com>

function addAggregated( filename, outfile)

	[f, t, v] = loadCsvFile( filename);
	ssz = rows(v);
	
	years = [];
	countryn = [];
	countries = [];
	for i=1:length(f)
		[num,status,strarray] = str2array(f{i}, ' ');
		years = [years num(1)];
		countryn = [countryn num(3)];
		countries = [countries {strarray{2}}];
	endfor
	
	#add aggregated by country
	for i=1:length(unique(countryn))
		pos = find ( countryn == unique(countryn)(i) );
		country = countries{pos(1)};
		count = sum( v( pos, :) );
		
		v = [v ; count];
		msg = sprintf("-1 %s %d", country, unique(countryn)(i) );
		f = [ f {msg} ];
	endfor
	
	#add aggregated by date
	for i=1:length(unique(years))
		pos = find ( years == unique(years)(i) );
		count = sum( v( pos, :) );
		
		v = [v ; count];
		msg = sprintf("%d <null> -1", unique(years)(i) );
		f = [ f {msg} ];
	endfor
	
	saveFile(f,t,v, outfile);
endfunction
