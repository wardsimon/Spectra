function mf_svfil
%
% MFIT function mf_loadr
% 		Get new output file
% 		MZ 29.11.94
%
hfile=findobj('tag','mf_OutFile');
hdir =findobj('tag','mf_OutDir');

%---------- Get current details on output file---------
outfile=get(hfile,'String');
outdir=get(hdir,'String');

%---- Dialogue to choose file -------------
[outfile, outdir]=uiputfile([outdir '*.mft'],'Select output file:');

if outfile~=0
	set(hfile,'String',outfile);
	set(hdir ,'String',outdir);
end;	

