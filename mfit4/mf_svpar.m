function  mf_svpar
% MFIT function mf_svpar
%	Save fitted parameters to file
% M. Zinkin 1.12.94
%----------------------------------------------------------------

%=============== Get information for output =====================

[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

if hmf_pars==0
	error('No parameter window open.')
end

%---------- Get pars, names, and errors --------------------------
h=get(hmf_pars, 'userdata');
pnames=[];
p=[];
std=[];
for i=1:size(h,1)
	p=[p; str2num(get(h(i,1),'string'))];
	s=max([0 str2num(get(h(i,2),'string'))]);
	std=[std; s];
	name=deblank(get(h(i,3),'string'));
	i=find(name==' ');
	name(i)='_'*ones(size(i));
	
	pnames=str2mat(pnames, setstr(name));
end
pnames(1,:)=[];

%----------- Time and Chi squared -------------------------------
t=fix(clock);
ChiSq=max([0 str2num(get(findobj('tag','ChiSq'),'String'))]);

%----------- Filenames and directories --------------------------
h=get(hmf_ctrl,'userdata');		

name=str2mat('DataFile', 'DataDir',...
				'FitFuncFile', 'FitFuncDir',...
				'OutFile', 'OutDir',...
				'LoadRoutineFile', 'LoadRoutineDir',...
				'FitRoutineFile', 'FitRoutineDir');

%============= Write to file  ====================================
s=['%%' setstr(abs('-')*ones(1,50)) '\n']; % make separator
outfile =[get(findobj('tag','mf_OutDir'),'String') ...
			 get(findobj('tag','mf_OutFile'),'String')];
fid=fopen(outfile,'a');
if fid==-1
	error('MFIT error: Bad output file');
end

%-------- Save header and variable values (eg. filenames) ---------
fprintf(fid,'%%MFIT Date %d.%d.%d  Time %d:%d:%d\n',t(3:-1:1),t(4:6));
fprintf(fid,s);
for i=1:size(name,1)
	h=findobj('tag',['mf_' deblank(name(i,:))]);
	par=get(h,'userdata');
	if isempty(par)
		par=get(h,'string');
	end
	fprintf(fid,'%s=%s\n',name(i,:),par);
end

%-------- save 'load' so that file can be run as a batch file ------
fprintf(fid,'load\n');

%-------- save parameter names and values --------------------------
fprintf(fid,s);
for i=1:length(p)
  fprintf(fid,'par %3d %s %14.6e  %13.6e\n',...
  				i,pnames(i,:),p(i),std(i));
end
fprintf(fid,s);

%-------- Save chi^2 ----------------------------------------------
fprintf(fid,'%% Norm. Chi^2 %12.3e\n%%\n',ChiSq);

fclose(fid);
    
