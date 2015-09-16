function mf_fcp(m)
%
% MATLAB file to change the fit control parameters
%
% DFM 24.3.96
%

%[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;
hmf_fcpdlg=findobj('Tag','mf_fcpdlg');
if nargin == 0, m = []; end


%========== Make dialog if it isn't there already =============

if isempty(hmf_fcpdlg)

%--------- Create figure window ----------------------------
   hmf_fcpdlg=figure('Position',[300 300 180 100],...
                        'MenuBar','none',...
      'Name','MFIT: Fit Control',...
      'NumberTitle','off',...
         'Color',get(0,'DefaultUicontrolBackgroundColor'), ...
      'Resize','off',...
      'Tag','mf_fcpdlg',...
      'Visible','off');

%--------- Make labels and edit boxes for limits------------

   boxh=20;
   NL = sprintf('\n');

   title=str2mat('dp','niter','stol');
   tag=str2mat('mf_dp','mf_niter','mf_stol');
   tip=str2mat(['Fit step, e.g. 0.l' NL '* Gradient: Partial derivative/increment step' NL '* Simplex: Exploration range in %'],'Max number of iterations, e.g. 20-100','Convergence tolerancy, e.g. 1e-5');
   pars = [ 0.005 20 1e-5];

   for i=1:3
      uicontrol(hmf_fcpdlg,...
            'Style','Text',...
            'String',title(i,:),...
            'ToolTipString', tip(i,:), ...
            'Position',[10 26+24*(3-i) 50 20]);
      uicontrol(hmf_fcpdlg,...
            'Style','Edit',...
            'BackgroundColor',[1 1 1],...
            'ForegroundColor',[0 0 0],...
            'Position',[60 26+24*(3-i) 100 20],...
            'String',pars(i),...
            'Tag',tag(i,:));
   end

%---------- Make ok and cancel buttons ---------------------

   uicontrol(hmf_fcpdlg,...
            'Style','push',...
            'String','Ok',...
            'Position',[25 2 50 20],...
            'Callback','mf_fcp(''set'')');
   uicontrol(hmf_fcpdlg,...
            'Style','push',...
            'String','Cancel',...
            'Position',[95 2 50 20],...
            'Callback','delete(gcf)');

   set(hmf_fcpdlg,'Visible','on');

%--------- else if exists bring to front ------------

elseif nargin == 0
   figure(hmf_fcpdlg)
end

if ~strcmp(m,'set')
   p=get(findobj('Tag','mf_fcpmenu'),'Userdata');
   h=[findobj('Tag','mf_dp   ');...
      findobj('Tag','mf_niter');...
      findobj('Tag','mf_stol ')];
   if length(h) == 3 & ~isempty(p)
    set(findobj('Tag','mf_dp   '),'String',num2str(p(1)));
    set(findobj('Tag','mf_niter'),'String',num2str(p(2)));
    set(findobj('Tag','mf_stol '),'String',num2str(p(3)));
   end
end

%--------- if self-call set new limits ---------------

if strcmp(m,'set')
   h=[findobj('Tag','mf_dp   ');...
      findobj('Tag','mf_niter');...
      findobj('Tag','mf_stol ')];
   for i=1:3
      p(i)=str2num(get(h(i),'String'));
   end
   set(findobj('Tag','mf_fcpmenu'),'Userdata',[p(1) p(2) p(3)]);
%   delete(gcf);
   if ~isempty(hmf_fcpdlg)
  set(hmf_fcpdlg,'visible','off');
   end
elseif ~isempty(findstr(m,'config'))
  j = findstr(m,'=');
  p = m((j(1)+1):length(m));
  p = str2num([ '[ ' p ' ]' ]);
  if ~isempty(p)
    h = findobj('Tag','mf_fcpmenu');
    if ~isempty(h)
      set(h,'Userdata',[p(1) p(2) p(3)]);
    end
    set(findobj('Tag','mf_dp   '),'String',num2str(p(1)));
    set(findobj('Tag','mf_niter'),'String',num2str(p(2)));
    set(findobj('Tag','mf_stol '),'String',num2str(p(3)));
  end
  set(hmf_fcpdlg,'Visible','off')
end
