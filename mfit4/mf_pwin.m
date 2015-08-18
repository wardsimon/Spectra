function [hmf_pars]=mf_pwin(name, pnames, p,redrawflag)
%
% MFIT function [hmf_pars]=mf_pwin(name, pnames, p)
%     Create new parameters window
%     MZ 29.11.94
%

if nargin < 3, p=[]; end
if nargin < 4, redrawflag=[]; end
funfile = [];

h = findobj('tag','mf_TextBoxHeight');
if isempty(h)
  Bheight = 18;
else
  Bheight =str2num(get(h,'string'));
end


Bspace=2;

[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;
npars=length(p);

if (isempty(p) | (p == 0))
  [npars,i] = size(pnames);
  p=zeros(npars,1);
  if npars == 0
    return
  end
end
if (hmf_pars~=0)
  oldnpars=size(get(hmf_pars,'userdata'),1);
  figure(hmf_pars);
  delete(gca);
else
  oldnpars=0;
end;

if isempty(redrawflag)
  redrawflag = 0;
end
if (nargin >3)
  oldnpars = 0;
end
if isempty(funfile)
  funfile = findobj('Tag','mf_FitFuncFile');
  if ~isempty(funfile)
    funfile = get(funfile,'string');
  else
    funfile='';
  end
end

if isempty(name) & ~isempty(funfile)
  name = funfile;
end
if ~isempty(name) & isempty(funfile)
  funfile = name;
end
if (length(name) > 15)
  name = [ name(1:13) '...' ];
end
if (isempty(hmf_pars) | hmf_pars==0 | npars~=oldnpars | redrawflag )        % No fig, or npars wrong

  if (hmf_pars~=0 & npars~=oldnpars)      % Fig exists, but npars wrong

    pos=get(hmf_pars,'position');
    set(hmf_pars,...
      'Visible','off',...
      'Name',[ 'MFIT: Parameters : ' name ],...
      'Position',[pos(1) pos(2) 323 (npars+2)*(Bheight+Bspace)]);
    delete(get(hmf_pars,'Children'));

  elseif isempty(hmf_pars) | hmf_pars==0          % No fig
    curfig=get(0,'CurrentFigure');
    h = findobj('tag','mf_ParWinPosition');
    if ~isempty(h)
      pos=sscanf(get(h,'string'),'%d');
    else
      pos = [ 30 70 ];
    end
    pos=[pos(1) pos(2) 323 (npars+2)*(Bheight+Bspace)];
    hmf_pars=figure( ...
      'tag','mf_ParWindow',...
      'HandleVisibility','on',...
      'Position',pos,...
      'Color',get(0,'DefaultUicontrolBackgroundColor'), ...
      'Name',[ 'MFIT: Parameters : ' name ],...
      'MenuBar','none',...
      'NumberTitle','off',...
      'Visible','off',...;
      'Resize','off');
    set(0,'CurrentFigure',curfig);
  end

  %----------- Create handles -------------------------------
  hpar=zeros(npars,1);
  hsig=zeros(npars,1);
  hpfix=zeros(npars,1);

  %-------- Make 'Parameter  Value  Error' titles: ----------------------
  title=str2mat('Parameter','  Value','Std Err');
  for i=1:3
      h = uicontrol(hmf_pars,...
              'Style','text',...
              'String',title(i,:),...
              'Position',[5+(i-1)*120 (npars+1)*(Bheight+Bspace)-Bspace 110 Bheight],...
              'ForegroundColor',[0 0 0]);
            if i == 1, set(h,'ToolTipString','Check parameters to fix their values'); end
  end

  for par=1:npars
    if isempty(pnames(par,:))
      pnames(par,:) = sprintf('par%i',par);
    end
   %------- Draw check boxes to fix/free parameters and parameter names-------
    pos=[ 20 (npars-par+1)*(Bheight+Bspace)+Bspace 115 Bheight];
    hpfix(par)=uicontrol(hmf_pars,...
                  'Style','checkbox',...
                  'ForegroundColor','k',...
                  'Position',pos,...
                  'HorizontalAlignment','left',...
                  'String',pnames(par,:));
    uicontrol(hmf_pars,...
                  'Style','text',...
                  'ForegroundColor','k',...
                  'Position',[ 0 (npars-par+1)*(Bheight+Bspace)+Bspace-2 16 Bheight ],...
      'string',sprintf('%i',par));
   % If p(par)==NaN box permanently checked. (Can't free parameter)
    if isnan(p(par))
    set(hpfix(par),...
      'Value',1,...
      'CallBack','set(gco,''Value'',1)');
    end

  %----------- Make editable text boxes for parameter values ------------------
    pos=[135 (npars-par+1)*(Bheight+Bspace)+Bspace 95 Bheight];
    hpar(par)=uicontrol(hmf_pars,...
                  'Style','edit',...
                  'BackgroundColor',[1 1 1],...
                  'ForegroundColor',[0 0 0],...
                  'String',num2str(p(par)),...
                  'HorizontalAlignment','right',...
                  'Position',pos);

  %-------- Draw text boxes for uncertainties in parameter values -----------
    pos=[234 (npars-par+1)*(Bheight+Bspace)+Bspace 85 Bheight];
    hsig(par)=uicontrol(hmf_pars,...
                  'Style','text',...
                  'BackgroundColor',[1 1 1],...
                  'ForegroundColor',[0.4 0.4 0.4],...
                  'HorizontalAlignment','right',...
                  'Position',pos);
  end

  % Draw fit results titles: 'Chi squared  Convergence'
  h=uicontrol(hmf_pars,...
        'Style','text',...
        'String','Chi^2',...
        'Position',[5 Bspace 50 Bheight],...
        'ForegroundColor',[0 0 0]);
  h=uicontrol(hmf_pars,...
        'Style','text',...
        'String','Q(Chi^2)',...
        'Position',[150 Bspace 70 Bheight],...
        'ForegroundColor',[0 0 0]);

  % Draw text boxes for fit results
  hfit(1)=uicontrol(hmf_pars,...
        'tag','ChiSq',...
        'Style','text',...
        'BackgroundColor',[1 1 1],...
        'ForegroundColor',[0.5 0.5 0.5],...
        'Position',[60 Bspace 75 Bheight],...
        'String','-----');
  hfit(2)=uicontrol(hmf_pars,...
        'tag','QChiSq',...
        'Style','text',...
        'BackgroundColor',[1 1 1],...
        'ForegroundColor',[0.5 0.5 0.5],...
        'Position',[223 Bspace 55 Bheight],...
        'String','-----');
else
  h=get(hmf_pars,'userdata');
  hpar=h(:,1);
  hsig=h(:,2);
  hpfix=h(:,3);
end

for par=1:npars
  set(hpfix(par),...
      'String',pnames(par,:),'ToolTipString', [ 'p(' num2str(par) ')=' pnames(par,:) ]);
end

set(hmf_pars,'userdata',[hpar hsig hpfix],...
    'Visible','on');

if ~isempty(name)
  set(hmf_pars,'name',[ 'MFIT: Parameters : ' name ]);
end
