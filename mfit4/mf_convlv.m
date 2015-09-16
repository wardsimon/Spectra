function [res] = mf_convlv(cmd, x, y)
% MFIT function for 1D convolution
% Convolution is done automatically as long as the window is opened
% and convolution check button is ON
% checks on dimension, center, norm, resampling  on X, etc is done
% automatically.
%
% Data for convolution may be:
% 1- undefined (at start)
% 2- transfered from the existing selected data set in MFit (name+X+Y)
% 3- transfered from an existing MView buffer (name+X+Y)
% 4- transfered the current fit function (name+parameters)
%
% The convolution data information is displayed in convlv window
%
% dependency: MFIT, funcs/convlv
% usage: mf_convlv(x,y)            use (x,y) for convlv
%        mf_convlv('data', x,y) use current Mfit data set (selection)
%        mf_convlv('mfit', 'data') use current Mfit data set (selection)
%        mf_convlv('mfit', 'func') use current Mfit fit function
%        mf_convlv('mview', NB)    use MView buffer #NB
%        mf_convlv('on')           activate convolution
%        mf_convlv('off')          un-activate convolution
%        mf_convlv('check')        perform checks and resampling
%        mf_convlv('create')       reate or show window
%        mf_convlv('hide')         hide window (cancel button). keep convlv data
%        mf_convlv('delete')       delete window and loose convlv data
%        mf_convlv('convlv', x, y) perform convolution of (x,y) by convlv data

% check of input parameters

% UserData structure

res = [];

if nargin == 0, cmd = []; end
[hmf_ctrl, hmf_data]=mf_figs;
if ~hmf_ctrl
  disp('mf_convlv: FATAL: requires MFit to run');
  return;
end

hmf_convlv = findobj('Tag', 'mf_convlv');

if isempty(cmd), cmd = 'convlv'; end
if findstr(cmd, 'convlv') & isempty(hmf_convlv), res = x; return; end
if isempty(hmf_convlv)
  cmd = [ cmd ' create' ];
elseif ~isempty(strfind(cmd, 'convlv'))
  UserData = get(hmf_convlv, 'UserData');
  if UserData.output.convlv & get(UserData.handles.checkbox, 'Value')
    res = convlv(x(:), UserData.output.y(:), 2);
  else res = x; end
  return
end % convlv

if ~ischar(cmd)
  if isnumeric(cmd) & isnumeric(x)
    mf_convlv('data',cmd, x); % transfert raw (x,y) data
    return;
  else
    disp('mf_convlv: warning: unknown non character argument');
    return;
  end
else
  cmd = lower(cmd);
end

NL = sprintf('\n');

if ~isempty(strfind(cmd, 'create')) | ~isempty(strfind(cmd, 'show'))
  if ~isempty(hmf_convlv)
    % raise/show
    set(hmf_convlv, 'Visible','on');
    UserData = get(hmf_convlv, 'UserData');
  else
    % create

    hmf_convlv = figure('Position',[800 200 300 200],...
       'Tag','mf_convlv',...
       'MenuBar','none',...
       'Name','MFit: 1D convolution',...
       'Resize','off',...
       'Visible','on',...
       'CloseRequestFcn','mf_convlv(''off hide'');');

   hmf_popup = uicontrol(hmf_convlv, ...
       'Style', 'popupmenu', ...
       'Callback','mf_convlv(''popup'');', ...
       'Tag','mf_convlv_popup', ...
       'ToolTipString',[ 'Define from here the convolution function' NL 'from MFit and MView' ], ...
       'String',{ 'Load...', ...
                  'MFit current data',...
                  'MFit current fit function', ...
                  'MView selected buffer'}, ...
       'Position',[ 5 175 100 20 ]);

   hmf_backgrd = uicontrol(hmf_convlv,...
       'Style','CheckBox',...
       'String','Bkg=0',...
       'Callback','mf_convlv(''background'');', ...
       'Tag','mf_convlv_background',...
       'Min',0, 'Max', 1, ...
       'Value', 1, ...
       'Position',[110 175 55 20], ...
       'ToolTipString',[ 'Check here to lower' NL ' the background to zero.' ]);

   hmf_center = uicontrol(hmf_convlv,...
       'Style','CheckBox',...
       'String','Center',...
       'Callback','mf_convlv(''center'');', ...
       'Tag','mf_convlv_center',...
       'Min',0, 'Max', 1, ...
       'Value', 1, ...
       'Position',[170 175 55 20], ...
       'ToolTipString',[ 'Check here to center' NL ' the convolution data' ]);

   hmf_norm = uicontrol(hmf_convlv,...
       'Style','CheckBox',...
       'String','Norm',...
       'Callback','mf_convlv(''normalize'');', ...
       'Tag','mf_convlv_norm',...
       'Min',0, 'Max', 1, ...
       'Value', 1, ...
       'Position',[230 175 50 20], ...
       'ToolTipString',[ 'Check here to normalize' NL ' the convolution data' ]);

   hmf_list = uicontrol(hmf_convlv, ...
       'Style', 'listbox', ...
       'Callback','mf_convlv(''check'');', ...
       'Tag','mf_convlv_listbox', ...
       'String',{ '[nothing loaded yet]', 'Select "Load From..." pop-up menu'}, ...
       'ToolTipString',[ '[nothing loaded yet]' NL 'Select "Load From..." pop-up menu' ], ...
       'Position',[ 5 30 280 140 ]);

    %------- dismiss and check box button -------------------------------------
   uicontrol(hmf_convlv,...
       'Style','PushButton',...
       'String','Dismiss',...
       'Tag','mf_convlv_dismiss',...
       'Position',[10 5 80 20], ...
       'CallBack', 'mf_convlv(''off hide'');', ...
       'ToolTipString',[ 'Click here to dismiss and' NL ...
                         'disable convolution' ]);

    uicontrol(hmf_convlv, ...
       'Style', 'pushbutton', ...
       'Callback','mf_convlv(''check plot'');', ...
       'ToolTipString','Click here to plot convolution data', ...
       'String','Show', ...
       'Position',[ 100 5 50 20 ]);

    hmf_checkbox = uicontrol(hmf_convlv,...
       'Style','CheckBox',...
       'String','Convolute',...
       'Callback','mf_convlv(''show toggle'');', ...
       'Tag','mf_convlv_convlv',...
       'Min',0, 'Max', 1, ...
       'Position',[160 5 110 20], ...
       'ToolTipString',[ 'Click here to enable' NL ' or disable convolution.' ]);

    UserData.input.type  = 'undefined';
    UserData.input.title = ''; % title to be displayed in mf_convlv window
    UserData.input.x     = []; % x or function parameters
    UserData.input.y     = [];
    UserData.input.name  = []; % file/source name or fit function name
    UserData.output.x    = []; % (x,y) data to be used for convolution
    UserData.output.y    = [];
    UserData.output.title= ''; % title to be displayed in mf_convlv window
    UserData.output.convlv= 0; % flag is true when convolution is possible

    UserData.output.center     = 1;
    UserData.output.background = 1; % flag is true when setting background to 0
    UserData.output.normalize  = 1;
    UserData.output.corrections= '';

    UserData.handles.listbox = hmf_list;
    UserData.handles.dialog  = hmf_convlv;
    UserData.handles.checkbox= hmf_checkbox;
    UserData.handles.popup   = hmf_popup;
    UserData.handles.center  = hmf_center;
    UserData.handles.background = hmf_backgrd;
    UserData.handles.normalize  = hmf_norm;

    set(hmf_convlv, 'UserData', UserData);

    cmd = [ cmd ' on' ];
  end
else
  UserData = get(hmf_convlv, 'UserData');
end % create/show

if ~isempty(strfind(cmd, 'popup'))
  value = get(UserData.handles.popup, 'Value');
  if value == 1, return; end
  switch value
  case 2
    mf_convlv('mfit','data');
  case 3
    mf_convlv('mfit','func');
  case 4
    mf_convlv('mview');
  end
  set(UserData.handles.popup, 'Value', 1);
  return
end

if ~isempty(strfind(cmd, 'hide'))
  set(hmf_convlv, 'Visible','off');
  return
end % hide

if ~isempty(strfind(cmd, 'delete'))
  delete(hmf_convlv);
  return
end % delete

if ~isempty(strfind(cmd, 'mfit'))
  if ~ischar(x)
    disp('mf_convlv: 2nd argument for MFit input data is not char');
    return;
  end
  data     = get(hmf_data,'userdata');
  datafile = get(findobj('Tag','mf_DataFile'),'String');
  datadir  = get(findobj('Tag','mf_DataDir'), 'String');
  if ~isempty(data)
    index=data(:,4);
    xs=data(:,1);
    ys=data(:,2);
  else index = []; xs = [];
  end
  if strcmp(x, 'data')
    if sum(index)< 5
      msg = sprintf('Error: MFit selected dataset does not contain enough points (%i) for convolution', sum(index));
      set(UserData.handles.listbox, 'ToolTipString', msg);
      mf_msg(msg);
      return
    end
    index = find(index & ~isnan(ys) & ~isinf(ys));
    UserData.input.x = xs(index);
    UserData.input.y = ys(index);
    UserData.input.name = [ datafile ' in ' datadir ];
    UserData.input.type = 'MFit data';
    UserData.input.title = '';
  else % end mfit data, begin mfit func
    fitfun=get(findobj('tag','mf_FitFuncFile'),'string');
    [p]=mf_rpars;
    if isempty(p)
      msg = 'Error: No fit function parameters defined yet in MFit';
      set(UserData.handles.listbox, 'ToolTipString', msg);
      mf_msg(msg);
      return;
    end
    [y, name, pnames]=feval(fitfun,[],p,1);
    UserData.input.x    = { fitfun, p };
    UserData.input.y    = { xs, ys };
    UserData.input.name = [ name ' (' fitfun ')' ];
    UserData.input.type = 'MFit function';

    to_add=cell(length(p),1);
    for index=1:length(p)
      to_add{index} = [ pnames(index,:) ' = ' num2str(p(index)) ];
    end
    UserData.input.title={ ...
      [ '1D Resolution function from ' UserData.input.type ], ...
      [ 'Source: ' UserData.input.name ], ...
      to_add{:}};
  end % mfit func
  x = [];
  cmd = [ cmd ' check' ];
end % mfit

if ~isempty(strfind(cmd, 'mview'))
  [hmv_ctrl, hmv_data]=mv_figs;
  if ~hmv_ctrl
    disp('mf_convlv: FATAL: requires MView to run');
    return;
  end
  buffers  =get(findobj('Tag','hmv_buffers'),'Userdata');
  tmv_radio=get(findobj('Tag','tmv_radio'),  'Userdata');
  if length(tmv_radio) ~= 1
    msg = 'Select ONE buffer to import in MView';
    set(UserData.handles.listbox, 'ToolTipString', msg);
    mf_msg(msg);
    return
  end
  oper_buffs=str2num(char(get(tmv_radio,'String')));
  buffer = buffers(oper_buffs);
  xi = buffer.xobs;
  yi = buffer.yobs;
  if ~length(xi) | ~length(yi)
     msg = sprintf('Selected MView buffer is invalid x(%i) y(%i)\n', ...
       length(xi), length(yi));
     set(UserData.handles.listbox, 'ToolTipString', msg);
     mf_msg(msg);
     return
  end
  UserData.input.x    = xi;
  UserData.input.y    = yi;
  UserData.input.name = buffer.datafile;
  UserData.input.type = [ 'MView buffer ' num2str(oper_buffs) ];
  UserData.input.title= '';
  cmd = [ cmd ' check' ];
end % mview

if ~isempty(strfind(cmd, 'data'))
  if ~length(x) | ~length(y) | length(x) ~= length(y)
    msg = sprintf('Raw input data is not consistent x(%i) y(%i)\n', ...
      length(x), length(y));
    set(UserData.handles.listbox, 'ToolTipString', msg);
    mf_msg(msg);
    return
  end
  UserData.input.x    = x;
  UserData.input.y    = y;
  UserData.input.name = '(x,y) data set';
  UserData.input.type = 'Raw data';
  UserData.input.title= '';
  cmd = [ cmd ' check' ];
end % data

if ~isempty(strfind(cmd, 'background'))
  h = UserData.handles.background;
  onoff = get(h, 'Value');
  UserData.output.background = onoff;
  cmd = [ cmd ' check' ];
end

if ~isempty(strfind(cmd, 'normalize'))
  h = UserData.handles.normalize;
  onoff = get(h, 'Value');
  UserData.output.normalize = onoff;
  cmd = [ cmd ' check' ];
end

if ~isempty(strfind(cmd, 'center'))
  h = UserData.handles.center;
  onoff = get(h, 'Value');
  UserData.output.center = onoff;
  cmd = [ cmd ' check' ];
end

if ~isempty(strfind(cmd, 'check')) | ~isempty(strfind(cmd, 'init'))
  UserData.output.convlv = 1;
  % set UserData.input.title based on stats
  if isempty(UserData.input.title) & isempty(strfind(UserData.input.type, 'function'))
    xi = UserData.input.x;
    yi = UserData.input.y;
    if ~sum(yi), return; end
    fmom=sum(xi.*yi)/sum(yi);
    smom=sqrt(sum(xi.^2.*yi)/sum(yi)-(sum(xi.*yi)/sum(yi))^2);
    UserData.input.title={ ...
      [ '1D Resolution function from ' UserData.input.type ], ...
      [ 'Source: ' UserData.input.name ], ...
      sprintf('Length %i points', length(xi)), ...
      sprintf('X=[%g:%g]', min(xi), max(xi)), ...
      sprintf('Y=[%g:%g]', min(yi), max(yi)), ...
      sprintf('Center=%g', fmom), ...
      sprintf('Width =%g (full)', smom), ....
    };
  end
  % get current MFit data and X step
  if nargin > 1
    mfit_x  = x;
  else
    mfit_x = [];
  end
  if isempty(mfit_x)
    data    = get(hmf_data,'userdata');
    if ~isempty(data)
      index   = data(:,4);    % selected data
      mfit_x  = data(:,1);
      index   = find(index);
      mfit_x  = mfit_x(index);
    end
  end
  if ~isempty(mfit_x);
    step_x  = abs(diff(mfit_x)); step_x  = step_x(find(step_x > 0));
  else
    mfit_x = []; step_x = [];
  end

  if ~length(step_x)
    msg = 'mf_convlv:MFit current data X sampling is void.';
    set(UserData.handles.listbox, 'ToolTipString', msg);
    fprintf(1, [ msg NL ]);
    fprintf(1, '          Can NOT convolute\n');
    UserData.output.convlv = 0;
    set(hmf_convlv, 'UserData', UserData);
    mfit_x = [];
  elseif ~length(mfit_x)
    fprintf(1, 'mf_convlv:MFit current data is empty.\n');
    fprintf(1, '          Can NOT convolute\n');
    UserData.output.convlv = 0;
    set(hmf_convlv, 'UserData', UserData);
  else
    step_x = mean(step_x);
  end
  % determine UserData.output

  if ~isempty(UserData.input.x) & ~isempty(mfit_x)
    % 1- first compute X axis for convolution data
    if ~isempty(strfind(UserData.input.type, 'function'))
      % first case: mfit func
      % xi is mfit x axis stored with the function
      xi = UserData.input.y{1};
      % yi = fitfun(xi, p)
      yi = feval(UserData.input.x{1}, xi, UserData.input.x{2});
    else
      % second case: data.
      % xi = min(input.x):mfit.step:max(input.x)
      xi = UserData.input.x;
      xi = min(xi):step_x:max(xi);
      % yi is interp(input.x, input.y, xi)
      yi = interp1(UserData.input.x, UserData.input.y, xi, 'cubic');
    end

    % 2- do checks
    if length(xi) < 5
      msg = sprintf('Convolution X data interpolated on current MFit X data has only %i points. step_x=%g\n', length(xi), step_x);
      set(UserData.handles.listbox, 'ToolTipString', msg);
      mf_msg(msg);
    end
    if length(xi) <= 3
      fprintf(1, '          Can NOT convolute\n');
      UserData.output.convlv = 0;
      set(hmf_convlv, 'UserData', UserData);
    end

    if UserData.output.convlv

      corrections = '';

      if UserData.output.center
        % 3- now center xi
        fmom    = sum(xi.*yi)/sum(yi);
        width_x = min( abs(fmom-min(xi)), abs(fmom-max(xi)) );
        width_x = max(floor(width_x/step_x), 2);
        width_x = width_x * step_x;

        xi = (fmom-width_x):step_x:(fmom+width_x);

        corrections = [ corrections ' centered' ];
      end

      % 4- then re-compute Y axis for convolution data and normalize yi
      if ~isempty(strfind(UserData.input.type, 'function'))
        % first case: mfit func
        % yi = fitfun(xi, p)
        yi = feval(UserData.input.x{1}, xi, UserData.input.x{2});
      else
        % second case: data.
        % yi is interp(input.x, input.y, xi)
        yi = interp1(UserData.input.x, UserData.input.y, xi, 'cubic');
      end

      if UserData.output.background
        yi = yi - min(yi);
        corrections = [ corrections ' backgrounded' ];
      end

      if UserData.output.normalize
        yi = yi./sum(yi);
        corrections = [ corrections ' normalized' ];
      end

      % 5- now compute basic stats and catenate info to listbox
      fmom=sum(xi.*yi)/sum(yi);
      smom=sqrt(sum(xi.^2.*yi)/sum(yi)-(sum(xi.*yi)/sum(yi))^2);
      UserData.output.title={ ...
        '', [ 'Convolution onto current MFit data:' ], ...
        sprintf('Length %i points', length(xi)), ...
        sprintf('X=[%g:%g]', min(xi), max(xi)), ...
        sprintf('Y=[%g:%g]', min(yi), max(yi)), ...
        sprintf('Center=%g', fmom), ...
        sprintf('Width =%g (full)', smom), ...
      };

      UserData.output.x = xi;
      UserData.output.y = yi;
      UserData.output.convlv = 1;
      UserData.output.corrections = corrections;

      % send title to listbox
      set(UserData.handles.listbox, 'String', { ...
        UserData.input.title{:}, ...
        UserData.output.title{:}, [ 'Corrections: ' corrections ] });
      set(UserData.handles.listbox, 'ToolTipString', ...
              [ '1D Resolution function from ' UserData.input.type NL ...
                'Source: ' UserData.input.name NL ...
                'Corrections: ' corrections ]);
    end % if UserData.output.convlv
  end % if ~isempty(UserData.input.x) & ~isempty(mfit_x)
  msg = 'Can not convolute (current convolution data is invalid).';
  set(hmf_convlv, 'UserData', UserData);
  if UserData.output.convlv, set(UserData.handles.listbox, 'ForegroundColor',[ 0 0 0 ]);
  else
    set(UserData.handles.listbox, 'ForegroundColor','red','ToolTipString',msg);
  end

  if UserData.output.convlv == 0 & get(UserData.handles.checkbox, 'Value')
    mf_msg(msg);
  end
end % check

if ~isempty(strfind(cmd, 'plot'))
  xi = []; xo = []; yi = []; yo = [];
  if ~isempty(UserData.input.x)
    if ~isempty(strfind(UserData.input.type, 'function'))
      xi = UserData.input.y{1};
      % yi = fitfun(xi, p)
      yi = feval(UserData.input.x{1}, xi, UserData.input.x{2});
    else
      xi = UserData.input.x;
      yi = UserData.input.y;
    end
    yi = yi /sum(yi);
  end
  if ~isempty(UserData.output.x)
    xo = UserData.output.x;
    yo = UserData.output.y;
  end
  if ~isempty(xi) | ~isempty(yo)
    f= figure('Name','MFit: 1D convolution: plot');
    hold on
    h = []; l = {};
    if ~isempty(xi)
      hi = plot(xi, yi, 'ro');
      l = { l{:}, 'Normalized Input data'};
    else hi = []; end
    if ~isempty(xo)
      ho = plot(xo, yo, 'bx');
      l = { l{:}, ['Convolution data' UserData.output.corrections]};
    else ho = [];  end
    h = [ hi ho ];
    legend(h, l);
    xlabel('MFit X axis'); ylabel('Signal');
    if isempty(xo) | UserData.output.convlv == 0
      title(['Convolution: 1D. Invalid data']);
    else
      title(['Convolution: 1D' UserData.output.corrections]);
    end
  end
end

if ~isempty(strfind(cmd, 'toggle'))
  hmf_convlv_convlv = UserData.handles.checkbox;
  onoff = get(hmf_convlv_convlv, 'Value');
  if gcbo == hmf_convlv_convlv, onoff = ~onoff; end
  if onoff, newstate = 'off' ;
  else      newstate = ' on'  ; end
  mf_convlv(newstate);
  return
end


if ~isempty(strfind(cmd, ' on'))
  hmf_convlv_convlv = UserData.handles.checkbox;
  set(hmf_convlv_convlv, 'Value', 1);
  h = findobj('tag','mf_mfit_convlv');
  if ~isempty(h)
    set(h,'Checked','on');
  end
end % on

if ~isempty(strfind(cmd, 'off'))
  hmf_convlv_convlv = UserData.handles.checkbox;
  set(hmf_convlv_convlv, 'Value', 0);
  h = findobj('tag','mf_mfit_convlv');
  if ~isempty(h)
    set(h,'Checked','off');
  end
end % off

if nargin > 0, res = UserData; end

set(hmf_convlv, 'UserData', UserData);
