function mf_opt(cmd)
%
% MFIT  function mf_opt(cmd)
%     Handle options: toggle grids and log/lin axes...
%     MZ 29.11.94 EF 04.97
%
% known commands : loglinx, logliny, drid, ebars,
%      autosave, autoguess, showext, about,
%      config, addtext, linksym, autochgx, revertfit
%      dobatch, rdpar
[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

if ~isempty(hmf_data) & hmf_data
  figure(hmf_data);
end

%---------- Toggle lin/log X-axis -----------------------
if strcmp(cmd,'loglinx')
   h=findobj('Tag','mf_loglinx');
   if strcmp(get(h,'Checked'),'off')
      set(gca,'Xscale','log')
      set(h,'Checked','on')
   else
      set(gca,'Xscale','lin')
      set(h,'Checked','off')
   end

%---------- Toggle lin/log Y-axis -----------------------
elseif strcmp(cmd,'logliny')
   h=findobj('Tag','mf_logliny');
   if strcmp(get(h,'Checked'),'off')
      set(gca,'Yscale','log')
      set(h,'Checked','on')
   else
      set(gca,'Yscale','lin')
      set(h,'Checked','off')
   end

%---------- Toggle grid on/off -----------------------
elseif strcmp(cmd,'grid')
   h=findobj('Tag','mf_gridonoff');
   if strcmp(get(h,'Checked'),'off')
      set(gca,'Xgrid','on')
      set(gca,'Ygrid','on')
      set(h,'Checked','on')
   else
      set(gca,'Xgrid','off')
      set(gca,'Ygrid','off')
      set(h,'Checked','off')
   end

%----------- Toggle error bar display on/off --------
elseif strcmp(cmd,'ebars')
   h=findobj('Tag','mf_ebarsonoff');
   hh=findobj('Tag','mf_ebars');
   if strcmp(get(h,'Checked'),'on')
     set(hh,'visible','off');
      set(h,'Checked','off');
   else
     set(hh,'visible','on');
      set(h,'Checked','on');
   end
   mf_uplot('err+sel');
elseif strcmp(cmd,'autoguess')
   h=findobj('Tag','mf_autoguess');
   hh=findobj('Tag','mf_AutoGuess');
   if strcmp(get(h,'Checked'),'on')
      set(hh,'string','0');
      set(h,'Checked','off');
   else
      set(hh,'string','1');
      set(h,'Checked','on');
   end
elseif strcmp(cmd,'autosave')
   h=findobj('Tag','mf_autosave');
   hh=findobj('Tag','mf_AutoSave');
   if strcmp(get(h,'Checked'),'on')
      set(hh,'string','0');
      set(h,'Checked','off');
   else
      set(hh,'string','1');
      set(h,'Checked','on');
   end
elseif strcmp(cmd,'showext')
   aload = get(findobj('Tag','mf_ExecAfterLoad'),'string');
   afit = get(findobj('Tag','mf_ExecAfterFit'),'string');
   prompt = {'After Load','After Fit','Any "mfit batch" command"'};
   def = {aload,afit,''};
   lineNo = 2;
   title = 'Enter optional  commands to execute';
   answer = inputdlg(prompt,title,lineNo,def);
   if ~isempty(answer)
    set(findobj('Tag','mf_ExecAfterLoad'),'string',answer{1});
    set(findobj('Tag','mf_ExecAfterFit'),'string',answer{2});
  if ~isempty(answer{3})
    mf_batch(answer{3});
  end
   end
elseif findstr(cmd,'About')
   helpstr = {cmd, ...
  'This program enables to',...
  '1- load data',...
  '2- view/zoom/select data',...
  '3- fit with any function or function set',...
  '4- choose the fit method',...
  '5- save your fit results and data',...
  '6- configure many features at users choice',...
  ' ',...
  'Authors : MZ, DFM, DAT, EF , ARW (1999)',...
  '*** MFit comes with ABSOLUTELY NO WARRANTY',...
  'This is free software, and you are welcome',...
  'to redistribute it under certain conditions',...
  'as specified in Licence files.'};
   helpdlg(helpstr,'MFit About');
elseif strcmp(cmd,'config')
   name = str2mat('AxesColor',...
        'FitColor',...
        'EbarColor',...
        'DataColor',...
        'DataColorSelected',...
        'MarkerSize',...
        'DataLineStyle',...
        'AxesFont', 'AxesFontSize',...
        'LabelFont', 'LabelFontSize',...
        'ExecAfterLoad',...
        'ExecAfterFit',...
        'IniFile',...
        'FitPoints');
   pars = [];
   for i=1:size(name,1)
  h=findobj('tag',['mf_' deblank(name(i,:))]);
  par=get(h,'userdata');
  if isempty(par) | size(par,1)>1
    par=get(h,'string');
  end
  if isempty(pars)
    pars = par;
  else
    pars = str2mat(pars,par);
  end
   end
   prompt = cellstr(name);
   def = cellstr(pars);
   lineNo = 1;
   title = 'MFit configuration variables';
   answer = inputdlg(prompt,title,lineNo,def);
   if ~isempty(answer)
  for i=1:size(prompt)
    h=findobj('tag',['mf_' char(prompt(i)) ]);
    set(h,'String',char(answer(i)));
  end
   end
elseif strcmp(cmd,'addtext')
   mf_msg('Click where to place new text...')
   refresh(hmf_data);
   figure(hmf_data);
   [tx,ty] = ginput(1);
% 'units','normalized',
   h=text(tx,ty,'',...
            'position',[0 1 0],...
            'Fontsize',10,...
            'FontName','courier',...
            'Tag','mf_text_user');

    h=mf_text(h,'');
    set(get(h,'userdata'),'visible','off');
    h = mf_text(h);
    set(h,'unit','data','position',[ tx ty 0]);

    mf_msg('done');
elseif strcmp(cmd,'linksym')
  aload = get(findobj('Tag','mf_ExecAfterLoad'),'string');
  if isempty(findstr(aload,'ssym'))
  mf_batch('ExecAfterLoad ssym(''fromfit+tomfit'');');
  end
  ssym('fromfit');
elseif strcmp(cmd,'autochgx')
   h=findobj('Tag','mf_autochgx');
   hh=findobj('Tag','mf_AutoRescale');
   if strcmp(get(h,'Checked'),'on')
      set(hh,'string','0');
      set(h,'Checked','off');
   else
      set(hh,'string','1');
      set(h,'Checked','on');
   end
elseif strcmp(cmd,'revertfit') & ~isempty(hmf_pars) & hmf_pars
   dpmf = get(findobj('Tag','ChiSq'),'Userdata');
   if ~isempty(dpmf)
   pin = dpmf(:,1);
   if size(dpmf,2) > 2
    pfix = dpmf(:,3);
    dp = dpmf(:,2);
   else
  pfix = dpmf(:,2);
  dp = [];
   end
   mf_upars(pin,dp,pfix);
   mf_gdata('noload');
   disp('Old parameters restored')
   else
  disp('No previous param set')
   end
elseif strcmp(cmd,'dobatch')
   str = get(findobj('tag','mf_CmdBatch'),'String');
   if ~isempty(str)
     mf_batch(str)
   end
elseif strcmp(cmd,'rdpar')
   % read par lines in Param section of MFT file
   [fullsections,filteredsections,sectionheaders] = mft_get('','Param','','par');
   if ~isempty(fullsections)
  mf_batch(fullsections);
   end
elseif strcmp(cmd,'viewout')
   outfile =get(findobj('tag','mf_OutFile'),'String');
   outdir  =get(findobj('tag','mf_OutDir'),'String');
   edit([outdir outfile]);
   mf_msg([ 'Edit output file ' outfile ]);
elseif strcmp(cmd,'viewdata')
   datafile =get(findobj('tag','mf_DataFile'),'String');
   datatdir =get(findobj('tag','mf_DataDir'),'String');
   edit([datatdir datafile]);
   mf_msg([ 'Edit data file ' datafile ]);
end

