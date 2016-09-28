function plotComparativeTempResiduAttrs(attr, files)

	if (length(files) == 0)
		error('Need to give me at least a filename to visualize');
	endif
	
	vs = [];
	[f,t,v1] = loadCsvFileMissing(files{1}, true);
	temp = v1(:, 1);	#extract temp information
	v1(:, 1) = [];
	vs = [ vs {v1}];
	temp = unique(temp);
	noyears = length(temp);
	noattrs = columns(v1);
	nocountries = rows(v1) / noyears;
	
	for pos=2:length(files)
		[f,t,v1] = loadCsvFileMissing(files{pos}, true);
		v1(:, 1) = [];
		vs = [ vs {v1}];
	endfor
	
	if (attr > 0 )
		i=1;
			legends = [];
			plott = [];
			for pos=1:length(files)
				legends = [ legends {sprintf('%s', files{pos})} ];
				plott = [ plott reshape(vs{pos}(:, attr), noyears, nocountries)(:,i) ];	
			endfor
		
			plot(temp, plott );
			title(sprintf('country: %s, attr %s', f{(i-1)*noyears + 1}, t{attr + 2}));
			ylabel(sprintf('attr'));
			xlabel("time");
			legend(legends);
	else
		for attr=1:noattrs
			i=1;
			legends = [];
			plott = [];
			for pos=1:length(files)
				legends = [ legends {sprintf('%s', files{pos})} ];
				plott = [ plott reshape(vs{pos}(:, attr), noyears, nocountries)(:,i) ];	
			endfor
		
			plot(temp, plott );
			title(sprintf('country: %s, attr %s', f{(i-1)*noyears + 1}, t{attr + 2}));
			ylabel(sprintf('attr'));
			xlabel("time");
			legend(legends)
			
			input(sprintf("Press any key for attr=%d\n", attr+1));
		endfor
	endif
	
endfunction
