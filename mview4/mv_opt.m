function mv_opt(cmd)
%
% MFIT  function mv_opt(cmd)
% 		Handle options: toggle grids and log/lin axes
% 		MZ 29.11.94
%
[hmv_ctrl, hmv_data]=mv_figs;

if hmv_data, figure(hmv_data); end

%---------- Toggle lin/log X-axis -----------------------
if strcmp(cmd,'loglinx')
   h=findobj('Tag','mv_loglinx');
   if strcmp(get(h,'Checked'),'off')
      set(gca,'Xscale','log')
      set(h,'Checked','on')
   else
      set(gca,'Xscale','lin')
      set(h,'Checked','off')
   end

%---------- Toggle lin/log Y-axis -----------------------
elseif strcmp(cmd,'logliny')
   h=findobj('Tag','mv_logliny');
   if strcmp(get(h,'Checked'),'off')
      set(gca,'Yscale','log')
      set(h,'Checked','on')
   else
      set(gca,'Yscale','lin')
      set(h,'Checked','off')
   end

%---------- Toggle grid on/off -----------------------
elseif strcmp(cmd,'grid')
   h=findobj('Tag','mv_gridonoff');
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
   h=findobj('Tag','mv_ebarsonoff');
   hh=findobj('Tag','mv_ebars');
   if strcmp(get(h,'Checked'),'on')
      set(hh,'visible','off');
      set(h,'Checked','off')
   else
      set(hh,'visible','on');
      set(h,'Checked','on')
   end

%----------- Fix axes limits or not? --------
elseif strcmp(cmd,'fixax')
   h=findobj('Tag','mv_fixax');
   if strcmp(get(h,'Checked'),'on')
      set(h,'Checked','off')
   else
      set(h,'Checked','on')
   end
   
elseif strcmp(cmd, 'legend')

  list = {'Best','None','BestOutside'};

  choice = listdlg('PromptString','Legend Style (Overlay)', ...
    'SelectionMode','single', ...
    'ListSize',[ 160 100 ],...
    'ListString',{'On (best)','Off (none)','Outside'});
  
 if ~isempty(choice)
   set(findobj('tag','mv_LegendStyle'),'string',list{choice});
 end
 mv_opt('legend_check');
 mv_msg('New settings will be applied for next Overlay');

elseif strcmp(cmd, 'legend_check')
  h = get(findobj('tag','mv_LegendStyle'),'string');
  if isempty(h), h = 'on'; end
  switch lower(h)
  case {'on','best'}
    h = 'on';
  case {'off','none'}
    h = 'none';
  case {'aside','outside','bestoutside'}
    h = 'outside';
  otherwise
    h = 'none';
  end
  set(findobj('Tag','mv_legend'),'Label',['Legend style: ' h]);
  set(findobj('tag','mv_LegendStyle'),'string',h);
elseif findstr(cmd,'About')
   helpstr = {cmd, ...
  'This program enables to',...
  '1- load data',...
  '2- view/zoom/select data',...
  '3- do basic operations between buffers',...
  '4- transfert data to MFit',...
  ' ',...
  'Authors : MZ, DFM, DAT, EF , ARW (1999)',...
  '*** MView comes with ABSOLUTELY NO WARRANTY',...
  'This is free software, and you are welcome',...
  'to redistribute it under certain conditions',...
  'as specified in Licence files.'};
   helpdlg(helpstr,'MView About');
  
end;


