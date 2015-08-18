function [i_start, i_end, s, pattern] = iGrep(s, pattern, option, engine)
% iGrep: regular expression search (ILL library)
%
% [i_start, i_end] = iGrep(string,pattern,option): regular expression search
% Searches for a regular expression 'pattern' in 'string'. If Java/mex is
% active, full regular expressions may be used. Else, only '*' and '?' wildcards
% are supported, and output is line numbers (option='l'). Then i_end is the 
% number of occurencies per line.
%
% string:  string to search in. May be a single char, a cellstr or char array
% pattern: pattern (reg. exp.)
% options: is an option string composed of any of:
%   o: searches the first occurency only (same as 'once' option)
%   c: output only the number of matches in s
%   i: ignore case
%   l: output the line number
%   f: file/path widlcards (e.g. /a/b/c/d*.ext). 
%       then replaces . by \., * by .* and ? by .? to obtain a valid RegExp.
%
% iGrep('help') opens the help web page about regular expression search
%
% See also: findstrs, strfind, findstr, strmatch
%
% Part of: iGrep utilities (ILL library)
% Author:  E. Farhi <farhi@ill.fr>. Feb 4th, 2003.

% Calls: gnu.regexp java class if available, grep_mex if available, regexp if available
% based on Java class gnu.regexp-1.1.4-dev (if Java is available)
% or reduced wildcard search (*, ?)

i_start = []; i_end = [];
if nargin == 1
  if strcmp(s, 'help') | strcmp(s, '--help') | strcmp(s, '-h')
    this_location = fileparts(which(mfilename));
    % calls GNU RegExp syntax help, with examples added
    this_location = [ this_location filesep 'doc' filesep 'REsyntax.html' ];
    if ~usejava('jvm') & ~exist('regexp') & (exist('grep_mex') ~= 3)
      this_location = [ this_location '#Simplified' ];
    end
    disp([ 'iGrep: URL=file:' this_location ]);
    web(this_location);
    return
  end
end
if nargin < 2
  disp('iGrep: requires at least 2 arguments (string, pattern)');
  return
end
if nargin < 3, option = ' '; end
if nargin < 4, engine = ''; end

flag.once = 0;
flag.once = 0;
flag.count= 0;
flag.icase= 0;
flag.lines= 0;
flag.file = 0;

if ~isempty(findstr(option, 'once')) options = 'o'; end
if ~isempty(find(option == 'o')) flag.once = 1; end
if ~isempty(find(option == 'c')) flag.count= 1; end
if ~isempty(find(option == 'i')) flag.icase= 1; end
if ~isempty(find(option == 'l')) flag.lines= 1; end
if ~isempty(find(option == 'f')) flag.file = 1; end
if ischar(s) & ~isempty(find(s == sprintf('\n')))
  s =strsplit(s);
  s = s(:);
end
if iscell(s) & ~iscellstr(s)
  index = find(~cellfun('isclass',s,'char'));
  s(index) = {' '};
end
s_in = s;
if flag.lines, 
  if iscellstr(s), s = char(s); end
  flag.lines=size(s, 2);
end
if iscellstr(s) | size(s,1) > 1, s = strjoin(s,'\n'); flag.lines= flag.lines+1;end
if ~ischar(s)
  disp('iGrep: acts on character strings/cellstr');
  return; 
end
if isempty(engine)
  if exist('regexp','builtin')   % matlab >= 6.5: use built-in regexp
    engine = 'regexp';
  elseif (exist('grep_mex') == 3) % use MeX regexp function
    engine = 'mex';
  elseif usejava('jvm') % Java enabled (matlab >= 6)
    engine = 'java';
  else
    engine = 'matlab';
  end % if else Java
end
if strcmp(engine, 'regexp') | strcmp(engine, 'java') | strcmp(engine, 'mex')
  if flag.file  % replaces . by \., * by .* and ? by .? to obtain a valid RegExp.
    pattern = strrep(pattern, '.','\.');  
    pattern = strrep(pattern, '*','.*');
    pattern = strrep(pattern, '?','.?');
  end
end
switch (engine)
case 'regexp'
  if flag.icase, f = 'regexpi'; else f = 'regexp'; end
  if flag.once, [i_start, i_end] = feval(f, s, pattern, 'once');
  else [i_start, i_end] = feval(f, s, pattern); end
case 'mex'
  [i_start, i_end] = grep_mex(s, pattern, option);
case 'java'
  % create the Regular Expression Java object
  RegExp = gnu.regexp.RE(pattern);
  if flag.icase, 
    RegExp = gnu.regexp.RE(pattern,RegExp.REG_ICASE+RegExp.RE_CHAR_CLASSES+RegExp.REG_MULTILINE); 
  end
  % only one match requested
  if flag.once, AllMatches= RegExp.getMatch(s);
  else AllMatches = RegExp.getAllMatches(s); end
  if ~flag.count
    if isempty(AllMatches), return; end
    i_start = zeros(1,length(AllMatches));
    i_end   = i_start;
    for index=1:length(AllMatches)
      ThisMatch = AllMatches(index);
      i_start(index) = ThisMatch.getStartIndex;
      i_end(index)   = ThisMatch.getEndIndex;
    end % for index
    % indexes start at 0 in Java
    i_start=i_start+1; i_end=i_end+1;
  end % if ~flag.count
  clear AllMatches RegExp
case 'perl'
  disp('iGrep: method perl is not supported')
otherwise % matlab method with reduced wildcards (* ?)
  engine = 'matlab';
  % no Java Virtual Machine available. Use reduced wildcard syntax.
  flag.lines = 1; flag.file=1; 
  if flag.file  % convert a real RE to a 'dir/ls' one
    pattern = strrep(pattern, '.*','*');
    pattern = strrep(pattern, '.?','?');
    pattern = strrep(pattern, '\.','.');
  end
  s = s_in;
  if flag.icase
    pattern = lower(pattern); s = lower(s);
  end
  s = deblank(cellstr(s));
  s = s(:);
  i_start= 1:length(s);
  i_end  = find(cellfun('isempty',s)); i_start(i_end) = 0;
  i_end  =i_start;
  p = pattern;
  while ~isempty(p)
    [t,p] = strtok(p,'*?'); % extract pattern parts
    if ~isempty(t)
      for index=1:length(s)
        if i_start(index) 
          s_tmp = s{index};
          tmp = strfind(s_tmp, t);
          if isempty(tmp), i_start(index) = 0; % remove non matching lines from the search
          else  % token is present
            if ~isempty(p)
              if (p(1) == '?')  % is the ? at the right place ?
                next_token = strtok(p, '*?');
                for i_tmp=1:length(tmp)
                  % next_token should follow 'tmp' position+'?' char
                  j_tmp = tmp(i_tmp)+length(t)+1;
                  if i_start(index) & ~isempty(next_token) & length(s_tmp) >= j_tmp
                    if ~strncmp(s_tmp(j_tmp:end), next_token, length(next_token))
                      tmp(i_tmp) = 0; % next token is absent
                    end
                  end
                end
                % all token occurencies where not followed by next token. 
                % remove line
                if isempty(find(tmp)), i_start(index) = 0; end
              end
            end
          end % if else token present
          if flag.lines, i_end(index) = length(tmp); end
        end % if valid index
      end % for
    end % ~isempty(t)
  end % while
  i_start = find(i_start); i_end=i_end(i_start);
  if flag.count, i_start = length(i_start);
  elseif flag.once & length(i_start)
    i_start = i_start(1); i_end=i_end(1);
  end
end % switch
if ~strcmp(engine, 'matlab')
  if flag.count, i_start = length(i_start);
  elseif flag.lines
    i_start = ceil(i_start/flag.lines);
    i_end   = ceil((i_end-1)/flag.lines);
  end
end
