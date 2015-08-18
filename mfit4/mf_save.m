function  mf_save(cmd)
% MFIT function mf_save
% Dialog to choose what to save
% can save : parameters, data+curve, ini (to mf_IniFile)
% other : getfile, go (selected before), 
%
% M. Zinkin 8.3.95
%
% Modified to allow fitted curve to be saved over a specified range of x with specified no. points
%   A.T. Boothroyd 4.07.02 (a.boothroyd1@physics.ox.ac.uk)

%----- If called from mfit (ie. not self-call) and the window isn't 
%		 already open, then make window.
		 
[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;
hmf_save = findobj('Tag','mf_savewhat');
outdir = pwd;
fitfun=get(findobj('tag','mf_FitFuncFile'),'string');
datafile = get(findobj('tag','mf_DataFile'),'string');

if (nargin==0 & isempty(hmf_save))

%====== Make dialog box if not there ===========================
	hmf_save=figure('Position',[200 200 200 160],...
				'Tag','mf_savewhat',...
				'Color',get(0,'DefaultUicontrolBackgroundColor'), ...
				'MenuBar','none',...
				'Name','MFIT: Save what:',...
				'NumberTitle','off',...
				'Resize','off',...
				'DefaultUicontrolForegroundColor',[0 0 0],...
				'Visible','off');

%---------- ok and cancel buttons --------------------
	uicontrol(hmf_save,...
            'Style','push',...
            'String','OK',...
            'Position',[30 2 60 18],...
            'Callback','mf_save(''go'')');
	uicontrol(hmf_save,...
            'Style','push',...
            'String','Cancel',...
            'Position',[110 2 60 18],...
            'Callback','set(gcf,''Visible'',''off'')');

%-------- Radio buttons --------------------------
	uicontrol(hmf_save,...							% Save fitted curve
				'Tag','mf_savecurve',...
				'Style','radio',...
				'Position',[30 95 160 20],...
				'String','Fitted curve',...
				'Value',0);
	uicontrol(hmf_save,...							% Save data
				'Tag','mf_savedata',...
				'Style','radio',...
				'Position',[30 115 160 20],...
				'String','Data',...
				'Value',0);
	uicontrol(hmf_save,...							% Save fit parameters and MFit main variables
				'Tag','mf_savefitpars',...
				'Style','radio',...
				'Position',[30 135 160 20],...
				'String','Fit Pars and Vars',...
				'Value',1);
   set(hmf_save,'Visible','on');
%--------- Text boxes to specify range and Npts for fit curve ------------------------
   h = findobj('tag','mf_TextBoxHeight');
    if isempty(h)
	    Bheight = 18;
    else
	    Bheight =str2num(get(h,'string'));
    end
    Bspace=2;
    xlimnames=str2mat('Npts','xmax','xmin');
    hax=get(hmf_data,'CurrentAxes');        % default is xmin, xmax for current plot and Npts = 100
    xlim=get(hax,'Xlim');  
    xlimvals=[100, max(xlim), min(xlim)];
    for i=1:3
	    uicontrol(hmf_save,...
	            'Style','text',...
	            'String',xlimnames(i,:),...
	            'Position',[40 35+(i-1)*(Bheight+Bspace) 60 Bheight],...
	            'ForegroundColor',[0 0 0]);
        pos=[90 35+(i-1)*(Bheight+Bspace) 60 Bheight];
        hh(i)=uicontrol(hmf_save,...
                    'tag',['mf_xlim',num2str(i)],...
	                'Style','edit',...
	                'BackgroundColor',[1 1 1],...
	                'ForegroundColor',[0 0 0],...
	                'String',num2str(xlimvals(i)),...
	                'HorizontalAlignment','right',...
	                'Position',pos);
    end
%===== Window exists so make current & visible ========================
elseif nargin==0
	set(findobj('Tag','mf_savewhat'),'Visible','on');
	figure(hmf_save);
	delete(gca);

%===== New output file name requested =================================
elseif strcmp(cmd,'getfile')

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

%===== Self-call - ok pressed on savewhat dialog ============================
elseif strcmp(cmd,'go')

%------ Work out what to save... --------------------------------------------
	svpars=get(findobj('Tag','mf_savefitpars'),'Value'); % fit pars button selected
	svcurv=get(findobj('Tag','mf_savecurve'),'Value');   % fit curve   "      "
	svdata=get(findobj('Tag','mf_savedata'),'Value');    % data        "      "

	if (isempty(svpars)) svpars = 1; end
	if (isempty(svcurv)) svcurv = 0; end
	if (isempty(svdata)) svdata = 0; end

%------ Save it...------------------------------------------------------------
	if  (svpars==1)
		mf_save('parameters')
	end

	if (svcurv==1 & svdata==1)
		mf_save('data+curve')
	elseif svdata==1
		mf_save('data')
	elseif svcurv==1
		mf_save('curve')
	end

%----- Turn window off -------------------------------------------------------
   set(findobj('Tag','mf_savewhat'),'Visible','off');


%===== Save parameters requested ===================================================

elseif strcmp(cmd,'parameters') & hmf_pars

%------ Get information for output 

   if isempty(hmf_pars)
   	errordlg('No parameter window open.','MFIT error:');
   	return
   end

%----- Get current parameters, names, and errors ------------------

   [p, dp, pfix,pnames]=mf_rpars;
   if isempty(dp)
	dp = 0*p;
   end

%------ Time and Chi squared 
   t=fix(clock);

%----- Filenames and directories 
   name=str2mat('DataFile', 'DataDir',...
   				 'FitFuncName','FitFuncFile', 'FitFuncDir',...
   				 'OutFile', 'OutDir');
   name=str2mat(name,...
   				 'LoadRoutineName','LoadRoutineFile','LoadRoutineDir',...
   				 'FitRoutineName','FitRoutineFile','FitRoutineDir');

   [hmf_ctrl, hmf_data, hmf_pars]=mf_figs;
   if ~isempty(hmf_data) & hmf_data
	hax = get(hmf_data,'CurrentAxes');
	hlx = get(findobj('Tag','mf_text_xlabel'),'String');
	if iscellstr(hlx), hlx = hlx{1}; end
	hly = get(findobj('Tag','mf_text_ylabel'),'String');
	if iscellstr(hly), hlx = hly{1}; end
	htl = get(findobj('Tag','mf_text_title'),'String');
	if iscellstr(htl), hlx = htl{1}; end
	[r2,rv,Chi] = mf_stats;
   else
	r2 = 0; rv = 0; Chi = [ 0 0 ];
	hax = []; hlx = []; hly = []; htl = []; hax = [];
   end
   
%--- Write to file...
   outfile = get(findobj('tag','mf_OutFile'),'String');
   outdir  = get(findobj('tag','mf_OutDir'),'String');
   if isempty(outdir)
	outdir = pwd;
   end
   fid=fopen([outdir filesep outfile],'a');
   if fid==-1
   	errordlg(['Bad output file ' outfile],'MFIT error');
   	return
   end

%-------- Save header and variable values (eg. filenames) 
   fprintf(fid,'%%MFIT Date %d.%d.%d  Time %d:%d:%d\n',t(3:-1:1),t(4:6));
   fprintf(fid,'%% Section : Vars - Data : %s - Function : %s .\n', datafile, fitfun);
   for i=1:size(name,1)
   	h=findobj('tag',['mf_' deblank(name(i,:))]);
   	par=get(h,'userdata');
   	if isempty(par) | size(par,1)>1
   		par=get(h,'string');
   	end
   	fprintf(fid,'%s = %s\n',name(i,:),par);
   end
   fprintf(fid,'%% [%s] is : [%s] versus [%s].\n',htl,hlx,hly);

%-------- save 'load' so that file can be run as a batch file 
   fprintf(fid,'load\n');

%-------- save parameter names and values 
   fprintf(fid,'%% Section : Parameters (%i) - Data : %s - Function : %s .\n', length(p), datafile, fitfun);
   for i=1:length(p)
     fprintf(fid,'par %3d %s %14.6e  %13.6e %3d\n',...
     				i,pnames(i,:),p(i),dp(i),pfix(i));
   end

%-------- Save stats : r2,rv,Chi
   fprintf(fid,'%% CorCoef %.3f -- RV %.3f -- ChiSq %.3f -- Q ChiSq %.3f\n', r2, rv, Chi(1), Chi(2));

   fclose(fid);
   mf_msg(['Saved to ' outdir filesep outfile]);
   disp(['Saved ' cmd ' to ' outdir filesep outfile]);


%===== Save data/fit requested ==========================================
elseif any([findstr(cmd,'data') findstr(cmd,'curve')]) & hmf_data

%------ Extract data from figure ----------------------
   hmf_data=findobj('tag','mf_DataWindow');
   if isempty(hmf_data)
   	mf_msg('No data window open.');
	return
   end

   data=get(hmf_data,'userdata');
   if isempty(data)
   	errordlg('No data available.','MFIT error:');
   end
   index=data(:,4);
   index = find(index ~= 0);
   data=data(index,1:3);
	coltitle=sprintf('%% %13s\t%13s\t%13s\t','x','y','err:y');
   mf_msg('Saving data...');

   %------ Evaluate the fit if necessary ---------------
   if ~isempty(findstr(cmd,'curve'))

   	fundir=get(findobj('tag','mf_FitFuncDir'),'string');

	hfit=findobj('Tag','mf_fitline');
	if ~isempty(hfit)
		yfit= get(hfit,'Ydata');
		xfit= get(hfit,'Xdata');
		if iscell(yfit)
			yfit = yfit{1};
			xfit = xfit{1};
		end
	end

   	if ~isempty(fitfun)
   		p=mf_rpars;
   		if ~isempty(p)
            if ~isempty(findstr(cmd,'data'))    % data + fit curve, so just evaluate
      				x=data(:,1);                % fit function at data points
				if (length(yfit) >= length(x))
					yfit = interp1(xfit,yfit,x);
				else
   	   	  			yfit = feval(fitfun,x,p);
				end
   		   		data=[data yfit];
   		   		coltitle=sprintf('%s%13s',coltitle,'yfit');
   			else                                % Save fit curve only (no data) over specified range and Npts.
                h1=findobj('tag','mf_xlim1');   % Npts
                h2=findobj('tag','mf_xlim2');   % xmax
                h3=findobj('tag','mf_xlim3');   % xmin
                x=linspace(str2num(get(h3,'string')),...
                           str2num(get(h2,'string')),...
                           str2num(get(h1,'string')))';
   	   	  		yfit = feval(fitfun,x,p);
               	data=[x yfit];
               	coltitle=sprintf('%13s\t%13s','xfit','yfit');
            end
   		end
   	end
   end

   %----------- Save the data ------------------------
   outfile = get(findobj('tag','mf_OutFile'),'String');	
   outdir  = get(findobj('tag','mf_OutDir'),'String');
   if isempty(outdir)
	outdir = pwd;
   end
   fid=fopen([outdir filesep outfile],'a');
   fprintf(fid,'%% Section : Data (%i)- Data : %s - Function : %s .\n', size(data,1),datafile, fitfun);
   fprintf(fid,'%s\n',coltitle);
   l=size(data,2);
   for i=1:size(data,1)
   	fprintf(fid,'%.6e\t',data(i,1:l-1));
   	fprintf(fid,'%.6e\n',data(i,l));
   end
   fprintf(fid,'%% End of Data %i lines, %i columns',size(data,1),l);
   fclose(fid);
   mf_msg(['Saved to ' outdir filesep outfile]);
   disp(['Saved ' cmd ' to ' outdir filesep outfile]);


%===== Save .ini requested ===================================================
elseif strcmp(cmd,'ini')
	
%----- Update control window position
   contpos=get(findobj('tag','mf_ControlWindow'),'position');
   set(findobj('tag','mf_ContWinPosition'),'string',sprintf('%d ',contpos(1:2)));

%----- Update parameter window position
   hpar=findobj('tag','mf_ParWindow');
   if ~isempty(hpar)
      parpos=get(hpar,'position');
      set(findobj('tag','mf_ParWinPosition'),'string',sprintf('%d ',parpos(1:2)));
   end

%----- Update figure window position and size
   hdata=findobj('tag','mf_DataWindow');
   if ~isempty(hdata)
      set(hdata,'unit','pixel');
      datpos=get(hdata,'position');
      set(findobj('tag','mf_FigurePosition'),'string',sprintf('%d ',datpos(1:2)));
      set(findobj('tag','mf_FigureSize'),'string',sprintf('%d ',datpos(3:4)));
   end

%----- Get .ini file name
	inifile=get(findobj('tag','mf_IniFile'),'string');
	if isempty(inifile), inifile = 'mfit.ini'; end
	fid=fopen(inifile,'w');
	if fid == -1
		fprintf(1,'%s could not be opened\n',inifile);
		return
	end

%----- Write {General} section
	disp('MFit variables');
	vars=get(findobj('tag','mf_VarList'),'userdata');     % variables to write
	fprintf(fid,'{General}\n');
	for i=1:size(vars,1)
		var=vars(i,:);
		var(isspace(var))=[];
		val=get(findobj('tag',['mf_' var]),'userdata');
		if isempty(val) | size(val,1)>1
			val=get(findobj('tag',['mf_' var]),'string');
		end
		fprintf(fid,'%s = %s\n',vars(i,:),val);
	end

%----- Write registry sections
   disp('MFit menus');
   SecTag = str2mat('Load Routines','Fit Functions','Fit Routines');
   NameTag = str2mat('mf_LoadRoutineName','mf_FitFuncName','mf_FitRoutineName');
   FileTag = str2mat('mf_LoadRoutineFile','mf_FitFuncFile','mf_FitRoutineFile');
   DirTag = str2mat('mf_LoadRoutineDir','mf_FitFuncDir','mf_FitRoutineDir');
   MenuTag = str2mat('mf_LoadRoutineMenu','mf_FitFuncMenu',  'mf_FitRoutineMenu');

   for i=1:size(SecTag(:,1))

      Name=get(findobj('tag',deblank(NameTag(i,:))),'userdata');
      Dir =get(findobj('tag',deblank(DirTag(i,:))),'userdata');
      File=get(findobj('tag',deblank(FileTag(i,:))),'userdata');
      hmenu = findobj(hmf_ctrl,'tag',deblank(MenuTag(i,:)));
      hmchild = get(hmenu,'children');
      
      fprintf(fid,'{%s}\n',deblank(SecTag(i,:)));

      for j=1:size(Name,1)
         if isempty(deblank(File(j,:)))
            fprintf(fid,'%s\n', Name(j,:));
         else
            fprintf(fid,'%s, %s, %s\n',File(j,:),Dir(j,:),Name(j,:));
         end
	 if j < size(Name,1) -1
           nextsep = get(hmchild(length(hmchild)-j-1),'separator');
	   if strcmp(nextsep,'on')
	     fprintf(fid,'==Separator==========\n');
	   end
	 end
      end
   end

%----- Write {Startup} section
	fprintf(fid,'{Startup}\n');

	if ~isempty(findobj('Tag','mf_mfunc'))
		disp('CFG:Multifunc');
		strmf = multifunc([],[],'getbatch');
		fprintf(fid,'%s',strmf);
	else
		strmf = [];
	end

	if ~isempty(findobj('Tag','hmf_wchx'))
		disp('CFG:X axis rescale');
		str = mf_chgx('getbatch');
		fprintf(fid,'%s',str);
	end

	if ~isempty(hmf_data) & hmf_data
		disp('CFG:load data');
    fprintf(fid,'load %s with %s\n', get(findobj('Tag','mf_DataFile'),'String'), get(findobj('Tag','mf_LoadRoutineFile'),'String'));
		fprintf(fid,'exec mf_gdata(''noload+nofit'');\n');
	end

	if strcmp(get(findobj('Tag','mf_AutoRescale'),'String'),'1')
		fprintf(fid,'exec mf_chgx(''newx'');\n')
	end

	if ~isempty(hmf_pars) & hmf_pars
		disp('CFG:Last parameter set')
		[p,dp,fixed] = mf_rpars;
		p = p(:)'; dp = dp(:)'; fixed = fixed(:)';
		p = num2str(p);
		dp = num2str(dp);
		fixed = num2str(fixed);
		str = sprintf('exec mf_upars([%s],[%s],[%s]);\n',p,dp,fixed);
		fprintf(fid,'%s',str);
	end
	hmf_fcpmenu =findobj('Tag','mf_fcpmenu');
	if ~isempty(hmf_fcpmenu)
		disp('CFG:Fit control parameters')
		p = get(hmf_fcpmenu,'Userdata');
		str = sprintf('exec mf_fcp(''config=[ %f %f %f ]'');\n',p);
		fprintf(fid,'%s',str);
	end

	if ~isempty(strmf)
		if findstr('spapp',strmf)
			if isempty(findstr(get(findobj('Tag','mf_ExecAfterLoad'),'String'),'ssym'))
				disp('CFG:MF ssym link')
				fprintf(fid,'exec ssym(''fromfit'');\n');
			end
		end
	end
   fclose(fid);
   fprintf(1,'%s/%s config file saved\n',pwd,inifile);
end


