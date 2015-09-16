function mf_gdata(cmd)
%
% MFIT function mf_gdata
%     Load new data file
%     MZ 29.11.94 EF 04.97
%

% cmds : nofit, noload, noplot
[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

%---------- Get load function details ---------------
loadfun=get(findobj('Tag','mf_LoadRoutineFile'),'string'); % load function file
datadir=get(findobj('Tag','mf_DataDir'),'string');           % data directory
datafile=get(findobj('Tag','mf_DataFile'),'string');         % data file
h = findobj('Tag','mf_PreviousDataFile');
if ~isempty(h)
  olddat = get(findobj('Tag','mf_PreviousDataFile'),'string');
else
  olddat = '';
end
if (nargin == 0)
  cmd = '';
end
if isempty(datadir), datadir = pwd; end
asv = str2num(get(findobj('tag','mf_AutoSave'),'string'));
if isempty(findstr(cmd,'noload')) & (asv == 1) & hmf_data
  if ~isempty(findstr(cmd,'direct')) & ~isempty(olddat)
    set(findobj('Tag','mf_DataFile'),'string',olddat);
  end
  mf_save('parameters');
  if ~isempty(findstr(cmd,'direct'))
    set(findobj('Tag','mf_DataFile'),'string',datafile);
  end
end
selected = 0;
if nargin==0 & isempty(findstr(loadfun,'frombase'))
  mf_msg([ 'Select file to load with ' loadfun ])
  datadir =get(findobj('Tag','mf_DataDir'),'String');
  dirsav = pwd;
  if ~isempty(datadir)
    if isdir(datadir), cd(datadir); end
  end
  
  % uigetfile on Windows can access the Net, but behaves strangly (prompt for opening/write)
  [datafile, datadir]=uigetfiles('*','Select data file','MultiSelect','off');

  cd(dirsav);
  if isempty(datafile) | ~ischar(datafile) datafile,
  return; end
  set(findobj('Tag','mf_DataDir'),'String',datadir);          % data directory
  set(findobj('Tag','mf_DataFile'),'String',datafile);
end

%datafile=get(findobj('Tag','mf_DataFile'),'string')
%datadir=get(findobj('Tag','mf_DataDir'),'string')

if ~isempty(hmf_data) & hmf_data
  userdata=get(hmf_data,'userdata');
  if ~isempty(userdata)
    x=userdata(:,1);
    y=userdata(:,2);
    err=userdata(:,3);
    selected=userdata(:,4);
  end

  xlab = get(findobj('tag','mf_text_xlabel'),'String');
  ylab = get(findobj('tag','mf_text_ylabel'),'String');
  titl = get(findobj('tag','mf_text_title'),'String');

  if isempty(xlab)
    xlab = get(get(findobj(hmf_data,'type','axes'),'Xlabel'),'string');
    if isempty(xlab)
      xlab = 'X';
    end
  end
  if isempty(ylab)
    ylab = get(get(findobj(hmf_data,'type','axes'),'Ylabel'),'string');
    if isempty(ylab)
      ylab = 'Y';
    end
  end

elseif ~isempty(findstr(cmd,'noload'))
%   disp('No user data');
    return;
end

if isempty(findstr(cmd,'noload'))
%------------ Call load function --------------------
  if ~isempty(findstr(cmd,'noprompt'))
    if ~isempty(hmf_data) & hmf_data
      userdata=get(hmf_data,'userdata');
      if ~isempty(userdata)
        selected = userdata(:,4);
      end
    end
    end

    mf_msg(['Loading ' datafile]);
    %loadfun
    %[datadir filesep datafile] = fullfile(datadir, datafile)
  [x y err xlab ylab]=feval(loadfun, fullfile(datadir, datafile));
  titl = [];

  filedata = dir(fullfile(datadir, datafile));
  if isempty(filedata), filedata = 0;  filedata.bytes=0; filedata.date=datestr(now); end
  set(findobj('Tag','mf_DataFile'),'ToolTipString', ...
    [ datafile ' from ' fullfile(datadir, datafile) sprintf('\n') ...
    filedata.date ' (' num2str(filedata.bytes) ' bytes) '  num2str(length(x)) ' points' sprintf('\n') ...
    '"' ylab '" vs. "' xlab '" ']);
end % load data

  x=x(:); y=y(:); err=err(:);

%------------ Exit if load unsuccessful -------------
  check=[size(x); size(y); size(err)];
  if any(any(~check))
    if size(err) == 0 & (size(x).*size(y) ~= 0)
      err = ones(size(y))*mean(y);
    else
      return;
    end
  end
if isempty(findstr(cmd,'noload'))
  mf_msg(['Loaded ' datafile]);
  fprintf(1,'* Loaded %s%s. %d points.\n', datadir, datafile, length(x));
  fprintf(1,'   [%s] vs [%s]\n',ylab,xlab);
end

%------------ Sort data ---------------------------
  i = find(~isnan(y) & ~isinf(y));
  x=x(i);
  y=y(i);
  if (length(err) ~= length(y))
    err = ones(size(y))*mean(y);
  else
    err=err(i);
  end

  [x, i]=sort(x);
  y=y(i);
  err=err(i);

%------------ Eliminate data with zero error -------
  i=find(err==0);
  if ~isempty(i)
    mf_msg('Data with zero error adjusted');
    if ~isempty(find(err))
      err(i)=min(err(find(err)));
    else
      err = ones(size(y))*mean(y)*100;
    end
  end
%------------ Make data window ---------------------

if ~isempty(findstr(cmd,'noload'))
  h=findobj('tag','mf_zoom_list');  % get old zooms
  oldzooms=get(h,'userdata');
  if iscell(oldzooms), oldzooms=[ min(x) max(x) min(y) max(y) ]; end
else
  oldzoom = [];
end

mf_msg('Plotting data');
if (~isempty(hmf_data) & hmf_data)
  figure(hmf_data);
  h=findobj(hmf_data,'type','text');
  delete(h)
  h=findobj(hmf_data,'type','line');
  delete(h)
end
if ~hmf_data
  hmf_data=mf_dwin(xlab, ylab);
end
figure(hmf_data);
%mf_text([],'mf_text_delete');

%---------- Attach data to userdata ------------------

[n,c]=size(x);
if c>n x=x'; end
[n,c]=size(y);
if c>n y=y'; end
[n,c]=size(err);
if c>n err=err'; end
if length(selected) ~= length(y)
  selected = ones(size(y));
end
if hmf_data
  set(hmf_data,'userdata',[x y err selected ]);
else
  disp('Can''t store data ! : no Data window')
end

%----------------- Set limits  --------------------
yrange=max(y)-min(y);
if yrange==0
   yrange=yrange+max([1e-3*mean(y) 1e-6]);
end
xrange=max(x)-min(x);
if xrange==0
   xrange=xrange+max([1e-3*mean(x) 1e-6]);
end

if isempty(findstr(cmd,'noload'))
  set(gca,'Xlim',[min(x)-0.01*xrange max(x)+0.01*xrange]);
  set(gca,'Ylim',[min(y)-0.02*yrange max(y)+0.02*yrange]);
end

%----------- Do the plot ------------------------------
%if isempty(findstr(cmd,'noplot'))
  hmf_data=mf_dwin(xlab, ylab);
%end

if ~isempty(titl)
  set(findobj('Tag','mf_text_title'),'String',titl);
end
h=findobj('Tag','mf_ebarsonoff'); % remove error bars if too many points
hh=findobj('Tag','mf_ebars');

if length(x) > 200 | mean(err) > mean(abs(y))/2
      set(hh,'visible','off');
      set(h,'Checked','off');
else
      set(hh,'visible','on');
      set(h,'Checked','on');
end


h = findobj('tag','mf_autochgx');
if ~isempty(h) & isempty(findstr(cmd,'noload'))
  if strcmp(get(h,'Checked'),'on')
    disp('Auto X axis rescale');
    mf_chgx('mfit+newx');
  end
end

h = findobj('tag','mf_ExecAfterLoad');
if ~isempty(h) & isempty(findstr(cmd,'noload'))
  todo = get(h,'String');
  if ~isempty(todo)
    fprintf(1,'Executing : %s\n',todo);
    eval(todo);
    if hmf_data
      set(hmf_data,'userdata',[x y err selected ]);
    else
      disp('Can''t store data ! : no Data window')
    end
  end
end
h = findobj('Tag','mf_PreviousDataFile');
if ~isempty(h)
  set(h,'string',datafile);
end

if isempty(findstr(cmd,'noplot'))
if ~isempty(findstr(cmd,'nofit'))
  mf_uplot('sel+err');
else
  mf_uplot('all');
end
end


if ~isempty(findstr(cmd,'noload')) & ~isempty(oldzooms)
  h=findobj('tag','mf_zoom_list');  % set old zooms
  set(h,'userdata',oldzooms);
end

set(hmf_data,'Name',[ 'MFIT:Data:' datafile ]);
