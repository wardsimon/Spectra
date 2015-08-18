function [data,datastr,scan]=multidata(file,matlabgui,choosedi)
% function [data,datastr,scan]=multidata(file,matlabgui,choosed)
%
% This function enables to load any data
% uses file selector and field chooser
% matlabgui = 1 for GUI window, 0 for none.
% choosed : set of field number to import, can be empty or not precised
%           - for main field (bigger field), * for all

% Author:  EF <manuf@ldv.univ-montp2.fr> 27.06.97
% Description:  General load routine

% uses : looktxt program

data=[];
choosed=[];
scan=[];
table=[];
datastr=[];
if nargin < 3, choosedi=[]; end
if nargin < 2, matlabgui=0; end

%===== Load file ============================================

%--------- Open data file ---------------------------------

[fid,m]=fopen(file,'r');
if (fid<0)
   disp([ 'File ' file ' could not be opened']);
   disp(m);
   return;
end

filestr = fread (fid);
filestr = setstr(filestr');
fclose(fid);

%------  looktxt analysis ---------------------------------

t=clock;
tmpfile = sprintf('lk%i',round(t(6)*10000));
cmd = [ '-F=' pwd filesep tmpfile ' -t -r -a -f -g ' file ];
if matlabgui cmd = [ cmd ' -v' ]; end
disp([ 'Starting file analysis : looktxt ' cmd ]);
%eval([ '! ' cmd ]);
looktxt(cmd);

nofile = 0;
verbose = 1;
if ~isempty(dir([ tmpfile '.m' ]))
  if verbose
    disp([ 'multiload: trying direct read in ' tmpfile ])
  end
  eval(tmpfile,'nofile = 1;');
end

if nofile
  [fid, message] = fopen([ tmpfile '.m' ], 'rt');

  if fid == -1
    disp(['multiload: error: looktxt could not analyse data file ' file ':' tmpfile ]);
    disp(['multiload: ' message ]);
    return
  else
    content = fread(fid, Inf);
    fclose(fid);
    if verbose
      disp([ 'multiload: using temporary file ' tmpfile '.m (' num2str(length(content)) ' bytes)' ]);
    end
    content = setstr(content');
    eval(content,[ 'disp(''multiload: could not evaluate data ' file ':' tmpfile '''); s=[]; ' ]);
  end
end

delete ([ tmpfile '.m' ])

% now table is 'table'
% [ id type start_pos end_pos rows columns ]
% with type = 1 (num) | 2 (text) | 64 (comment)
if isempty(table)
  disp('No data was found in file')
  return
end
test = bitand(1,table(:,2));
numtable = table(find(test),:);
[n,c] = sort(numtable(:,5).*numtable(:,6));
numtable = numtable(c,:);
[n,c] = size(numtable);
numtable = numtable((n:-1:1),:); % revert order : maxid is first
[maxel,maxid] = max(numtable(:,5).*numtable(:,6));
if (n > 30)
  fprintf(1,'Warn : Many (%i) fields... \n',n);
  fprintf(1,'Press Ctrl-C to abort if needed\n\n');
end

endheaderpos = min([ numtable(maxid,3) numtable(maxid,4) ]);
if (endheaderpos > 500)
  endheaderpos = min([ numtable(1,3) numtable(1,4) ]);
end

header = filestr(max(1,endheaderpos-500):endheaderpos);
fprintf(1,'Header : ...\n%s\n',header);

if (n>1)
  fprintf(1,'Found %i fields\n',n);

  if (matlabgui & isempty(choosedi))

  liststring = [];

  for i=1:n
        toadd = sprintf('n%i (array %ix%i)',numtable(i,1),numtable(i,5),numtable(i,6));
  if (i == maxid)
    toadd = [ toadd ' [Biggest]' ];
  end
  liststring = strvcat(liststring, toadd);
  end
  [choosed,okx] = listdlg('ListString',cellstr(liststring),...
  'PromptString','Select parts to import','InitialValue',maxid,...
  'Name','Enter field to import','ListSize',[300 160]);
  if ~okx
    choosed = [];
  end

%============ Now wait until ok or cancel pressed =========================


  else % no matlabgui or choosed in input args
  for i=1:n
    fprintf(1,'%i - n%i = (%ix%i)\n',i,numtable(i,1),numtable(i,5),numtable(i,6));
  end
  if (isempty(choosedi) | strcmp(choosedi,'-')) choosed=maxid;
  elseif strcmp(choosedi,'*') choosed = 1:n;
  else choosed = choosedi; end
  if isstr(choosed) choosed = str2num(choosed); end
  end % end if matlabgui
else
  fprintf(1,'Found one field... ok\n');
  choosed = 1;
end

if (isempty(choosed))
  disp('Nothing imported.');
  return;
end

data = [];
n = length(choosed);
fprintf(1,'Getting %i fields : ',n);

nr = sum(numtable(choosed,5));
nc = max(numtable(choosed,6));
data = zeros(nr,nc);
r=1;
for i=1:n
  nr = numtable(choosed(i),5);
  nc = numtable(choosed(i),6);
  ni = ['n' num2str(numtable(choosed(i),1))];
  fprintf(1,'%s  ',ni);
  eval([ 'nf = ' ni ';' ]);
  [nr,nc] = size(nf);
  data((r:(r+nr-1)),1:nc) = nf;
  r=r+nr;
end
[nr,ncol] = size(data);
fprintf(1,'\nData matrix has %ix%i elements.\n',nr,ncol);

datastr='';
choosed=choosed(:)';
scan=num2str(choosed);
