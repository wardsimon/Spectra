function rdir = getfiles(filename, interactive)
% getfiles: returns a list of file names to import.
%    files = getfiles(filename, interactive)
%    returns a cellstr of full file names to import.
%    it handles wildchars '*' and directory names
%    input:  filename=single char name, or char array of names, or cellstr
%              array. URLs and ZIP/GZIP are supported (a local uncompressed copy
%              of the file is created)
%            interactive=0 (default)|1 for asking user choice
%    output: files   =single name as char, or a cell string array of names
%
% See also: uigetfiles, geturl
%
% Part of: iFiles utilities (ILL library)
% Author:  E. Farhi <farhi@ill.fr>. Feb 4th, 2003.

% calls:   iGrep, iFiles/geturl, iFiles/uigetfiles

rdir = {};
if nargin < 1, return; end
if nargin < 2, interactive = []; end
if isempty(interactive), interactive = 0; end

if ~ischar(filename) & ~iscellstr(filename)
  disp('iFiles/getfiles: "filename" should be a char or cell string array.')
  return
end

if isempty(filename)
  if interactive
    [filename, filepath] = uigetfiles('', 'LOAD: Select a data file');
      if ~isempty(filename)
        filename = strcat(filepath,filesep,filename);
      else
        return  % user cancelled choice
      end
  else return
  end % if interactive
end % if no file name

if iscellstr(filename)
  if length(filename) == 1
    rdir = char(filename);  % single cell
  else  % cell str array
    for index = 1:length(filename)
      tmp = getfiles(filename{index}, interactive);
      if iscellstr(tmp), rdir = { rdir{:}, tmp{:} };
      else rdir = { rdir{:}, tmp }; end
    end
    rdir = rdir(:); % as a clomun
    return
  end
end % if iscellstr

if ischar(filename) & size(filename,1) > 1
  % char array
  rdir = getfiles(cellstr(filename), interactive);
  return
end

% from there filename is a single char, and is not empty
[filepath,name,ext,versn]=fileparts(filename);
% if compressed file, force to call Java method
if ~isempty(ext) & ~isempty(strmatch(lower(ext),{'.gz','.gzip','.tgz','.zip'}))
  if ~strncmp(filename,'file:',5) & ~strncmp(filename,'http://',7) ...
                                  & ~strncmp(filename,'ftp://',6)
    if isempty(filepath), filename = [ pwd filesep filename ]; end
    filename = [ 'file:' filename ];
  end
else  % local uncompressed file: use normal procedure
  if strncmp(filename,'file:',5), filename = filename(6:end); end
end
% URL and GZIP/ZIP handling
if strncmp(filename,'http://',7) | strncmp(filename,'ftp://',6) ...
                                 | strncmp(filename,'file:',5)
  filename = geturl(filename);
  if isempty(filename), return; end
end

if isdir(filename)
  filename = [ filename filesep '*' ]; % get full directory content
end
[filepath,name,ext,versn]=fileparts(filename);  % 'file' to search

if ~isempty(find(filename == '*'))  % wildchar !!
  if isempty(filepath), filepath = pwd; end
  this_dir = dir([ filepath filesep filename]);
  if isempty(this_dir), this_dir = dir(filepath); end
  index = find(real(char(this_dir.isdir)) == 0);
  this_dir = char(this_dir.name);
  this_dir = (this_dir(index,:));
  rdir = cellstr(this_dir); % original full directory listing as cell
  % we search 'regex with *' in these files
  r = [ name ext versn ];
  index = iGrep(rdir, r, 'fl');  % filenames, and want line numbers
  % add directory path to file names
  rdir = strcat([ filepath filesep ], char(rdir{index}));
  if interactive & size(rdir,1)>1
    [select, ok] = listdlg('ListString', cellstr(rdir),...
      'SelectionMode','multiple','ListSize',[ 300 160 ], ...
      'Name',[ 'Import files ' filename], 'PromptString',...
      {'Select the file(s) to load from wildcard:',filename });
    if ok == 0
      rdir = {};
      return
    else
      rdir = rdir(select,:);
    end
  end
  rdir = cellstr(rdir);
else
  rdir = filename;
end

if iscellstr(rdir) & length(rdir) == 1
  rdir = char(rdir);
end