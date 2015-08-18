function  mf_svdat(cmd)
% MFIT function mf_svdat
% Function to save data and fitted curve to ascii file
% Possible values of cmd are:
%  
%
% M. Zinkin 8.3.95

hmf_data=findobj('tag','mf_DataWindow');

%------ Extract data from figure ----------------------
if isempty(hmf_data)
	error('MFIT error: No data window open.');
end

data=get(hmf_data,'userdata');
if isempty(data)
	error('MFIT error: No data available.');
end
index=data(:,4);
index = find(index ~= 0);
data=data(index,1:3);

mf_msg('Saving...');

%------ Evaluate the fit if necessary ---------------
if ~isempty(findstr(cmd,'curve'))

	fitfun=get(findobj('tag','mf_FitFuncFile'),'userdata');
	fundir=get(findobj('tag','mf_FitFuncDir'),'userdata');

	if fitfun~=''
		p=mf_rpars;
		if ~isempty(p)
         if ~isempty(findstr(cmd,'data'))       % data+curve, so just evaluate
   			x=data(:,1);                        % fit function at data points
	   	  	yfit=feval(fitfun,x,p);
		   	data=[data yfit];
		   else                                   % just curve, so generate 100
		      hax=get(hmf_data,'CurrentAxes');    % points between current axis
            xlim=get(hax,'Xlim');               % limits
            x=linspace(min(xlim),max(xlim),100)';
            data=[x feval(fitfun,x,p)];
         end
		end
	end

end

%----------- Save the data ------------------------
outfile =get(findobj('tag','mf_OutFile'),'String');	
outdir  =get(findobj('tag','mf_OutDir'),'String');	
fid=fopen([outdir outfile],'a');
l=size(data,2);
for i=1:size(data,1)
	fprintf(fid,'%.6e\t',data(i,1:l-1));
	fprintf(fid,'%.6e\n',data(i,l));
end
fclose(fid);

mf_msg(['Saved to ' outdir outfile]);

