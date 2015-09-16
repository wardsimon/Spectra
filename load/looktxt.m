% [main, root, file, content] = looktxt('filename [options]') Import text data
% Converts any text file into a Matlab text file
% Then, type 'filename' to import data.
% 'main'    is the main numeric field number.
% 'root'    is the variable name.
% 'file'    is the file name.
% 'content' is either the full data as a structure,
%              or the string of the output file for direct evaluation
%
% Use "looktxt" or "looktxt -h" for more help about options.
% Example :
%          looktxt foo.txt -p="." -s="\t\v,;" -c="#"
% E.Farhi 01/97 and K. Dabertran 04/99. v0.3 (06/11/00)

% Author:  EF <farhi@ill.fr>
% Description: Import any text file data. Need 'looktxt.c' source code.

% this code is used when no MeX is found
function [mainfield, rootname, file, content] = looktxt(args)
if nargin == 0
        args = '';
end
mainfield = -1;
rootname = [];
file = [];
content = [];
table = [];


  t=clock;
  tmpfile = sprintf('lk%i',round(t(6)*10000));
  if isempty(args), args = ' '; end
  args = [ args ' -t -r -F=' tmpfile ];

if exist('matlabroot')
  system([ 'looktxt  ' args ],'');
  tmpfilem = [ tmpfile '.m' ];
  [fid, msg] = fopen(tmpfilem,'r');
  if (fid ~= -1)
    F = fread(fid, Inf);
    fclose(fid);
    content = setstr(F');
    rootname = tmpfilem;
    delete(tmpfilem);
    file = tmpfile;
    mainfield = 0;
  end
elseif exist('OCTAVE_VERSION')
  if isempty(findstr(args, '-o=oct')) & isempty(findstr(args, '-o="oct"')) & isempty(findstr(args, '-S'))
    args = [ args ' -S' ];
  end
  [s,w] = system([ 'looktxt ' args ]);
  disp(s);
  if ~isempty(findstr(args, '-o=oct')) | ~isempty(findstr(args, '-o="oct"'))
    tmpfilem = [ tmpfile '.oct' ];
    load('-force',tmpfilem);
    content = [ 'load -force ' tmpfilem ];
    rootname = tmpfilem;
  else
    if ~isempty(findstr(args, '-S'))
      eval(tmpfile,'s=[];');
      if ~isempty(s)
        content = s;
        table = s.table;
      end
    end
    tmpfilem = [ tmpfile '.m' ];
    rootname = tmpfilem;
  end
  file = tmpfile;
  if ~isempty(dir(tmpfilem)), system([ 'rm ' tmpfilem ]); end
end

if ~isempty(table) & (mainfield == -1)
  index = find(table(:,2) == 1);
  mn = table(index,5).*table(index,6);
  [mnmax, mnmaxi] = max(mn);
  mainfield = table(index(mnmaxi),1);
end
