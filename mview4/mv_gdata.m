function mv_gdata(cmd)
%
% MVIEW function to get data into mview
%
% DFM 26.2.95
%
% Modified to take monitor as a variable as well
% ARW 22.11.00

%----- Unload Names of files

loadfun=get(findobj('Tag','mv_LoadRoutineFile'),'String');
loaddir=get(findobj('Tag','mv_LoadRoutineDir'),'String');


if nargin==0
  mv_msg('Select file(s) to load')

  if (isunix), mask = '*';
  else mask = '*.*'; end
  datadir =get(findobj('Tag','mv_DataDir'),'String');
  dirsav = pwd;
  if ~isempty(datadir)
    if isdir(datadir), cd(datadir); end
  end
  %[datafile, datadir]=uigetfiles(mask, 'MVIEW: Select data file(s)');
  [datafile, datadir]=uigetfiles(mask, 'MVIEW: Select data file(s)','MultiSelect','on');
  cd(dirsav);
  if isempty(datafile), return; end
  if isreal(datafile)
    if datafile(1) == 0, return; end
  end
  
  filename = fullfile(char(datadir),char(datafile));
elseif ischar(cmd)
  filename = cmd;
elseif isreal(cmd)
  datadir =get(findobj('Tag','mv_DataDir'),'String');
  datafile=get(findobj('Tag','mv_DataFile'),'String');
  filename = fullfile(datadir,datafile);
end

filename = getfiles(filename); % handle RegExp and check dirs

if iscellstr(filename)
  for index=1:length(filename)
    this_file = filename{index};
    if index > 1, this_file = [ this_file, ',file_list' ]; end
    mv_gdata(this_file);
  end
  return
elseif ischar(filename)
   if ~isempty(findstr(filename, ',file_list'))
    is_serie = 1;
    filename = strrep(filename, ',file_list', ''); % removes the list flag
   else is_serie = 0; end
   filename = getfiles(filename);
   [datadir, name, ext, ver] = fileparts(filename);
   datafile = [ name, ext, ver];
   set(findobj('Tag','mv_DataDir'),'String',datadir);    % data directory
   set(findobj('Tag','mv_DataFile'),'String',datafile);
   if is_serie, datafile1 = [ datafile ',file_list' ];
   else datafile1 = datafile; end
else
  disp('mv_gdata: unknown type for file name');
  return
end

%------------ Call load function --------------------

mv_msg(['Loading ' datafile]);
dir_old=cd;
%cd(loaddir);
x=[];y=[];err=[];mon=[];
eval('[x y err x_label y_label mon]=feval([loadfun], fullfile(datadir,datafile1));',...
     '[x y err x_label y_label]=feval([loadfun], fullfile(datadir,datafile1));');
%cd(dir_old);

%------------ Exit if load unsuccessful -------------

check=[size(x); size(y); size(err)];
if any(any(~check))
   return;
end
mv_msg(['Loaded ' datafile]);

g_label=datafile;

%----- Check if monitor was passed, if not create it.
if isempty(mon), mon=ones(size(x)); end;

%----- Sort into increasing x values

[x,perm]=sort(x);
y=y(perm);
err=err(perm);
mon=mon(perm);

%----- Plot data

mv_graph(x,y,err,x_label,y_label,g_label)

%----- Load data into buffers

incr_buffs=mv_ubuff(x,y,err,mon,x_label,y_label,g_label,datafile);
set(findobj('Tag','hmv_CurrentBuffer'),'Userdata',incr_buffs);

%----- Change text in window to file name
filedata = dir(fullfile(datadir,datafile));
hmv_text=get(findobj('Tag','hmv_text'),'Userdata');
set(hmv_text(incr_buffs),'String',g_label);
if ~isempty(filedata), 
   set(hmv_text(incr_buffs),'ToolTipString',[ datafile ' from ' fullfile(datadir,datafile) sprintf('\n') ...
     filedata.date ' (' num2str(filedata.bytes) ' bytes) ' num2str(length(x)) ' points' sprintf('\n') ...
     '"' y_label '" vs. "' x_label '"' ]);
end

return

