function mv_xylim(m)
%
% MATLAB function to set limits for a graph in MVIEW
%
%
[hmv_ctrl, hmv_data]=mv_figs;
hax=get(hmv_data,'CurrentAxes');
hmv_lims=findobj('Tag','mv_xylimdlg');

%========== Make dialog if it isn't there already =============

if (hmv_data~=0 & isempty(hmv_lims) & nargin==0)

%--------- Create figure window ----------------------------

   hmv_lims=figure('Position',[300 300 180 130],...
                   'MenuBar','none',...
                   'Name','MFIT: Axis limits',...
                   'NumberTitle','off',...
                   'Color',get(0,'DefaultUicontrolBackgroundColor'), ...
                   'Resize','off',...
                   'Tag','mv_xylimdlg',...
                   'Visible','off');

%--------- Make labels and edit boxes for limits------------

   boxh=20;
   lims=[get(hax,'Xlim') get(hax,'Ylim')];
   title=str2mat('X min','X max','Y min','Y max');
   tag=str2mat('mv_xlim','mv_xlim','mv_ylim','mv_ylim');
   
   for i=1:4
      uicontrol(hmv_lims,...
               'Style','Text',...
               'String',title(i,:),...
               'Position',[10 26+24*(4-i) 50 20]);
      uicontrol(hmv_lims,...
               'Style','Edit',...
               'String',num2str(lims(i),6),...
               'BackgroundColor',[1 1 1],...
               'ForegroundColor',[0 0 0],...
               'Position',[60 26+24*(4-i) 100 20],...
               'Tag',tag(i,:));
   end

%---------- Make ok and cancel buttons ---------------------

   uicontrol(hmv_lims,...
            'Style','push',...
            'String','Ok',...
            'Position',[25 2 50 20],...
            'Callback','mv_xylim(''set'')');
   uicontrol(hmv_lims,...
            'Style','push',...
            'String','Cancel',...
            'Position',[95 2 50 20],...
            'Callback','delete(gcf)');

   set(hmv_lims,'Visible','on');

%--------- else if exists bring to front ------------

elseif (~isempty(hmv_lims) & nargin==0)

   figure(hmv_lims)

%========= if self-call set new limits ===============

elseif strcmp(m,'set')

   h=[findobj('Tag','mv_xlim'); findobj('Tag','mv_ylim')];
   for i=1:4
      p(i)=str2num(get(h(i),'String'));
   end
   set(hax,'Xlim',[min(p(1:2)) max(p(1:2))])  
   set(hax,'Ylim',[min(p(3:4)) max(p(3:4))])
   delete(gcf);

else

   m=strrep(m,',',' ');
   p=sscanf(m,'%f');
   if length(p)>1
      set(hax,'Xlim',[min(p(1:2)) max(p(1:2))])  
   end
   if length(p)>3
      set(hax,'Ylim',[min(p(3:4)) max(p(3:4))])  
   end

end
