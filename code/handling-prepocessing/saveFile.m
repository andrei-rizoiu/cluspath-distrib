function saveFile( lineHeader, columnHeader, vals, file )

	% function that will save in a data file in the text format.
	% note: line and column headers can be absent
	
	[fid, MSG] = fopen(file, 'w');
	
	if ( fid == -1)
		MSG
		return;
	endif

	% write the column headers, if they exist
	if (columns(columnHeader) > 0 )
		if ( columns(columnHeader) == columns(vals) )
			fprintf(fid, "\t");
		endif
		for i=1:columns(columnHeader)-1
			fprintf(fid, "%s\t", columnHeader{i});
		endfor
		fprintf(fid, "%s\n", columnHeader{columns(columnHeader)});
	endif
	
	%write each line
	for i=1:rows(vals)
		if ( columns(lineHeader) > 0 )
			fprintf(fid, "%s\t", lineHeader{i});
		endif
		fprintf(fid, "%d\t", vals(i, 1:(end-1) ));
		fprintf(fid, "%d\n", vals(i, end ));
	endfor
	
	fclose(fid);
end
