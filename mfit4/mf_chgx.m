function [newx, rescaletable, xpol,x,y,peak] = mf_chgx(x,y,peak,cmd)
% [newx, rescaletable, pol ,x,y,peak] = mf_chgx(x,y,peak,cmd)
%
% Change X axis with reference points.
% cmd:
%     hide, show, newx, getpol, ginput, max, peak, transfert
%     getbatch
%     xpol=<polynome coefficients>

h_data = [];
if nargin == 1 & isstr(x)
  cmd = x;
  x = [];
end
name = '';
newx = []; rescaletable = []; pol = [];
if nargin < 4, cmd = []; end
if nargin < 3, peak= []; end
if nargin < 2, y= []; end
if nargin < 1, x= []; end
if isempty(cmd)
  cmd = 'show';
end
hsp_ctrl = findobj('tag','hsp_ctrl');
hmf_data = findobj('Tag','mf_DataWindow');
toactivate = get(0,'CurrentFigure');
if ~isempty(findstr(cmd,'mfit'))
  h_data = hmf_data;
  if isempty(peak)
    peak = get(findobj('Tag','mf_autoguess'),'userdata');
  end
  name = 'MFit';
elseif ~isempty(findstr(cmd,'ssym'))
  hsp_plot = []; hsp_app = [];
  hsp_ctrl = findobj('tag','hsp_ctrl');
  hsp_pars = findobj('tag','hsp_pars');
  userdata = get(hsp_ctrl,'userdata');
  if ~isempty(userdata)
    h_data  = userdata(1);
  end
  peak=str2num(get(findobj(hsp_pars,'tag','hsp_ppeak'),'String'));
  toactivate = hsp_ctrl;
  name = 'SSym';
else
  h_data = findobj('Tag','hmf_wchx');
  peak = get(findobj('Tag','hsp_xpk'),'userdata');
  name = get(findobj('Tag','hsp_xmod'),'String');
end

if ~strcmp(get(findobj('Tag','hsp_xmod'),'String'),name)  % mode has changed
  cmd = [ 'force+' cmd ];
end

if isempty(name)
  name = get(findobj('Tag','hsp_xmod'),'String');
end

if strcmp(name,'MFit')
  toactivate = hmf_data;
elseif strcmp(name,'SSym')
  toactivate = hsp_ctrl;
end

% params input

if isempty(x)
  userdata = get(h_data,'userdata');
  if ~isempty(userdata)
    x = userdata(:,1);
  else
    x = [];
  end
end

if isempty(y)
  userdata = get(h_data,'userdata');
  if ~isempty(userdata)
    y = userdata(:,2);
  else
    y = [];
  end
end

if isempty(peak) & ~isempty(y) & ~isempty(x)
  disp('Compute peak matrix')
  peak = mf_fndpks(x,y);
  [dummy, peak] = mf_rmspks(x,y,peak);
end

% build window if necessary

if ~isempty(findstr(cmd,'force'))
  delete(findobj('tag','hmf_wchx'));
end

if isempty(findobj('tag','hmf_wchx'))

  Fposl=200; Fposb=400; Fwidth=300; Fheight=160;
  hmf_wchx = figure('Position',[Fposl Fposb Fwidth Fheight],...
    'Color',get(0,'DefaultUicontrolBackgroundColor'), ...
    'MenuBar','none',...
    'Name',[ name ' Change Axis' ],...
    'NumberTitle','off',...
    'Resize','off',...
    'userdata',[],...
    'Tag','hmf_wchx', ...
    'Visible','off');

  NL = sprintf('\n');

  uicontrol(hmf_wchx,...    % Incoming
    'Tag','hsp_xinc',...
    'Position',[15 50 60 90 ],...
    'Style','edit',...
    'Min',0,'Max',2,...
    'String','',...
    'FontName','Times','FontSize',12,...
    'ForegroundColor',[0 0 0],...
    'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','right');

  uicontrol(hmf_wchx,...    % OldX
    'Tag','hsp_xold',...
    'Position',[150 50 60 90 ],...
    'Style','edit',...
    'Min',0,'Max',2,...
    'String','',...
    'FontName','Times','FontSize',12,...
    'ForegroundColor',[0 0 0],...
    'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','right');

  uicontrol(hmf_wchx,...    % NewX
    'Tag','hsp_xnew',...
    'Position',[ 220 50 60 90 ],...
    'Style','edit',...
    'Min',0,'Max',2,...
    'String','',...
    'FontName','Times','FontSize',12,...
    'ForegroundColor',[0 0 0],...
    'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','right');

% text ------------------------------------------------

  uicontrol(hmf_wchx,...    % Incoming
    'String','Incoming',...
    'FontName','Times','FontSize',12,...
    'Position',[ 15 145 60 15 ],...
    'Style','text',...
    'ToolTipString',[ 'Some X axis values from the original data set' NL 'May be refined with nearest maximum ands peak' ],...
    'ForegroundColor',[1 1 1],...
    'HorizontalAlignment','left');

  uicontrol(hmf_wchx,...    % Old X
    'String','Old Axis',...
    'FontName','Times','FontSize',12,...
    'Position',[ 150 145 60 15 ],...
    'Style','text',...
    'ToolTipString','X axis values from the original axis',...
    'ForegroundColor',[1 1 1],...
    'HorizontalAlignment','left');

  uicontrol(hmf_wchx,...    % New X
    'String','New Axis',...
    'FontName','Times','FontSize',12,...
    'Position',[ 220 145 70 15 ],...
    'Style','text',...
    'ToolTipString','X axis values to replace original ones',...
    'ForegroundColor',[1 1 1],...
    'HorizontalAlignment','left');

% buttons ------------------------------------------------

  uicontrol(hmf_wchx,...    % get Points
    'Style','PushButton',...
    'String','Get Points',...
    'ToolTipString','Get X axis values from data window (mouse click)',...
    'Callback','mf_chgx(''ginput'');',...
    'Position',[10 4 70 20]);
    uicontrol(hmf_wchx,...    % Max
    'Style','PushButton',...
    'String','Max',...
    'ToolTipString','Refine Incoming Values to the nearest maxima position',...
    'CallBack','mf_chgx(''max'');', ...
    'Position',[ 5 30 40 20]);
    uicontrol(hmf_wchx,...    % Peak
    'Style','PushButton',...
    'String','Peak',...
    'ToolTipString','Refine Incoming Values to the nearest peak position',...
    'CallBack','mf_chgx(''peak'');', ...
    'Position',[47 30 40 20]);
    uicontrol(hmf_wchx,...    % Transfert
    'Style','PushButton',...
    'String','Transfert =>',...
    'ToolTipString','Tranfert Incoming values into the Original ones for axis conversion',...
    'CallBack','mf_chgx(''transfert'');', ...
    'Position',[80 70 70 20]);

  uicontrol(hmf_wchx,...    % OK
    'Style','PushButton',...
    'String','Ok',...
    'Callback','mf_chgx(''hide+newx'');',...
    'Position',[ 193 4 50 20]);
    uicontrol(hmf_wchx,...    % cancel
    'Style','PushButton',...
    'String','Cancel',...
    'CallBack','mf_chgx(''hide'');', ...
    'Position',[ 245 4 50 20]);

  uicontrol(hmf_wchx,...    % Order
    'Style','popup',....
    'Tag','hsp_xord',...
    'String','Select|Order1|Order 2|Order 3|Order 4',...
    'Value',1,...
    'ToolTipString','Interpolation polynome order',...
    'Callback','mf_chgx(''getpol'');',...
    'Position',[ 235 30 60 20 ]);

  uicontrol(hmf_wchx,...    % polynome
    'Tag','hsp_xpol',...
    'Position',[ 100 30 130 22],...
    'Style','edit',...
    'String','',...
    'ToolTipString','Interpolation polynome',...
    'FontName','Times','FontSize',12,...
    'ForegroundColor',[0 0 0],...
    'BackgroundColor',[1 1 1]);

% invisible -------------------------------------

  uicontrol(hmf_wchx,...
    'Tag','hsp_xpk',...
    'userdata',[],...
    'Visible','off');

  uicontrol(hmf_wchx,...
    'Tag','hsp_xmod',...
    'userdata',[],...
    'Visible','off');


else
  hmf_wchx = findobj('Tag','hmf_wchx');
end

% commands ------------------------------------------------

if findstr(cmd,'hide')
  set(hmf_wchx,'Visible','off');
end

if findstr(cmd,'show')
  set(hmf_wchx,'Visible','on');
  figure(hmf_wchx);
end

% now get data
xinc = str2num(get(findobj('Tag','hsp_xinc'),'String'));
xold = str2num(get(findobj('Tag','hsp_xold'),'String'));
xnew = str2num(get(findobj('Tag','hsp_xnew'),'String'));
xpol = str2num([ '[ ' get(findobj('Tag','hsp_xpol'),'String') ']' ]);
xord = get(findobj('Tag','hsp_xord'),'Value') -1;

if ~isempty(findstr(cmd,'xpol'))
  j = findstr(cmd,'=');
  xpol= cmd((j(1)+1):length(cmd));
  xpol = str2num([ '[ ' xpol ' ]' ]);
end

if ~isempty(findstr(cmd,'ginput'))
  figure(toactivate);
  drawnow
  fprintf(1,'Activating %s\n',get(toactivate,'Name'));
  disp('Hit a Key or Button 2-3 to end Point Collecting')
  but = 1;
  nb = 0;
  fprintf(1,'X: ');
  while (but == 1)
    drawnow
    [xin,yin,but] = ginput(1);
    if (but == 1)
      xinc = [ xinc(:) ; xin ];
      nb = nb + 1;
      fprintf(1,'%.3f ',xin);
    end
  end
  fprintf(1,'- %i points\n',nb);
  figure(hmf_wchx);
end

if ~isempty(findstr(cmd,'max'))
  for i = 1:size(xinc,1)
    [dummy,imax] = min(abs(x(peak(:,1)) - xinc(i)));
    if ~isempty(imax)
      xinc(i) = x(peak(imax,1));
    end
  end
  xinc = unique(xinc);
end

if ~isempty(findstr(cmd,'peak'))
  for i = 1:size(xinc,1)
    [dummy,imax] = min(abs(peak(:,4) - xinc(i)));
    if ~isempty(imax)
      xinc(i) = peak(imax,4);
    end
  end
  xinc = unique(xinc);
end

if ~isempty(findstr(cmd,'transfert'))
  xold = [ xold(:) ; xinc(:) ];
  xinc = [];
end

if ~isempty(findstr(cmd,'getpol'))
  if isempty(xnew)
    if ~isempty(xold) & ~isempty(xpol)
      xnew = polyval(xpol,xold);
    else
      disp('Can''t interpolate : enter "New Axis" values first...')
    end
  else
    if xord == 0
      if length(xpol) > 1
        xord = length(xpol) -1;
      else
        xord = 2;
      end
    end
    if length(xnew) <= xord
      disp('Not enough points for interpolation')
    else
      fprintf(1,'Interpolation order %i\n',xord);
      xpol = polyfit(xold,xnew,xord);
      fprintf(1,'P = [ '); fprintf(1,'%.3g ',xpol); fprintf(1,' ]\n');
    end
  end
end

if ~isempty(findstr(cmd,'newx'))
  if length(xpol) > 1 & size(peak,1) > 1
    newx = polyval(xpol,x);
    peak(:,4) = polyval(xpol,peak(:,4));
    % should also transfert to main app : mfit, ssym
    if ~isempty(newx) & ~all(newx == mean(newx))
      if ~isempty(hmf_data) & hmf_data
        tomfit(newx,y);
        mf_rbbox('zoomreset');
      end
      if ~isempty(hsp_ctrl) & hsp_ctrl
        ssym(newx,y);
      end
      x = newx;
    end
  else
    disp('Polynome isn''t defined. Can''t change axis.')
  end
end

set(findobj('Tag','hsp_xinc'),'String',num2str(xinc,6));
set(findobj('Tag','hsp_xold'),'String',num2str(xold,6));
set(findobj('Tag','hsp_xnew'),'String',num2str(xnew,6));
set(findobj('Tag','hsp_xpol'),'String',num2str(xpol,'%.8g '));
set(findobj('Tag','hsp_xord'),'Value',1);
set(findobj('Tag','hmf_wchx'),'userdata',[ x(:) y(:) ]);
set(findobj('Tag','hsp_xpk'),'userdata',peak);
set(findobj('Tag','hsp_xmod'),'String',name);

if ~isempty(findstr(cmd,'getbatch'))
  if ~isempty(xpol)
    newx = num2str(xpol);
  else
    newx = '1 0';
  end
  newx = sprintf('exec mf_chgx(''mfit+xpol=%s'');\n',newx);
  return
end

if nargout == 0
  newx = [ cmd ' done' ];
end

