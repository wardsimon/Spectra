function [fullsection,sectstrmat,fullsecstr] = mft_get(file,sectiontype,sectionnumber,sectionfilter)
% [fullsections,filteredsections,sectionheaders] = mft_get(file or string,sectiontype,filter1,filter2)
% reads an .mft file (or any section delimited file/string)
% search for 'Section' and sectiontype in file,
% then select those containing filter1 (can also be 'all' or ask user if empty)
% and eventually look for filter2 lines.
% exemple :
%     [..] = mft_get('toto.mft','Param','all',['Temp','Intensity']);
%   will extract lines containing Temp and Intensity in all sections
%   identified by the word 'Param'.
% See also : mf_gpar, much simpler, in order to extract 2 parameters

% EF 05/98

% we collect in file all section headers.
% then search users choice
% and available parameters names
if nargin < 4, sectionfilter=''; end
if nargin < 3, sectionnumber=''; end
if nargin < 2, sectionfilter=''; end
if isempty(sectiontype)
  sectiontype = 'Parameters';
end

fullsection=[];
sectstrmat = [];

eol = [ 10 13 ]; % end of line characters

% get file -------------------------------
if isempty(file)
  [file, datadir]=uigetfiles('*.mft','Select data file','MultiSelect','off');
  if file==0 return; end
end

[fid, msg] = fopen(file,'r');
if fid == -1
  s = file;
else
  F=fread(fid);
  s=setstr(F');
  fclose(fid);
end

sectindexes = findstr(s,'Section');
if isempty(sectindexes)
  sectindexes = findstr(s,sectiontype);
end

if isempty(sectindexes)
  if fid ~= -1
    fprintf(1,'Sorry : no Section or %s found in file %s\n',sectiontype,file);
  else
    fprintf(1,'Sorry : no Section or %s found in input string\n',sectiontype);
  end
  return
end

% reject contiguous sections beginings

sectindexesdiff = [ 2 diff(sectindexes) ];
j = find(sectindexesdiff > 1);
sectindexes = sectindexes(j);

fullsecstr = '';
fullsecidx = [];

for i= sectindexes
  if ~isempty(sectiontype)
    j = findstr(s(i:(i+30)),sectiontype);
  else
    j = 1;
  end
  if ~isempty(j)
    j = i+j(1)-1;
    sechdr = strtok(s(i:length(s)),eol); % this is the section header line.
    fullsecstr = strvcat(fullsecstr,sechdr);
    fullsecidx = [ fullsecidx j ];
  end
end

% now look if some sections are contiguous, then gather.

fullsecidxdiff = [ 2 diff(fullsecidx) ];
j = find(fullsecidxdiff > 1);

for i=find(fullsecidxdiff < 2);
  k = find(j < i);
  if ~isempty(k)
    k = k(length(k));
    fullsecidx(i) = k;
  else
    fullsecidx(i) = 0;
  end
end

if isempty(fullsecidx)
  disp('No section found')
  return;
end

if isstr(sectionnumber) & findstr(sectionnumber,'all')
  sectionnumber = 1:length(fullsecidx);
end

if iscellstr(sectionnumber)
  sectionnumber = char(sectionnumber);
end

if isstr(sectionnumber)
  outsecnumber = [];
  for i=1:length(fullsecidx)
    for j=1:size(sectionnumber,1)
      if fullsecidx(i) & findstr(fullsecstr(i,:),sectionnumber(j,:))
        outsecnumber = [ outsecnumber i ];
      end
    end
  end
  sectionnumber = outsecnumber;
end


if isempty(sectionnumber) & size(fullsecstr,1) == 1
  sectionnumber = 1;
end

if isempty(sectionnumber)
  [sectionnumber,ok] = listdlg('ListString',cellstr(fullsecstr),'PromptString','Select sections to import','Name','File scanning','ListSize',[300 160]);
  if ~ok
    return
  end
end

i = find(sectionnumber > 0 & sectionnumber <= length(fullsecidx));
sectionnumber = sectionnumber(i);
fullsection = '';
sectionstart = fullsecidx(sectionnumber);
fullsecstr = fullsecstr(sectionnumber,:);
fullsecidx = [ fullsecidx length(s) ];
sectstrmat = '';

for i = 1:length(sectionnumber)
  sstart = sectionstart(i);
  send = find(sectindexes > sstart);
  if isempty(send)
    send = length(s);
  else
    send = sectindexes(send(1))-1;
  end
  if sstart
    fullsection = sprintf('%s%% Section : %s\n',fullsection,s(sstart:send));

% now convert into strmat to search filter

    tmpsect = s(sstart:send);
    [line, tmpsect] = strtok(tmpsect,eol);
    sectstrmat = strvcat(sectstrmat,[ '% Section : ' fullsecstr(i,:) ]);

    while (~isempty(line) & ~isempty(tmpsect))
      line = line(1:length(line-1));
      if isempty(sectionfilter)
        sectstrmat = strvcat(sectstrmat,line);
      else
        for j=1:size(sectionfilter,1)
          if ~isempty(findstr(line,sectionfilter(j,:)))
            sectstrmat = strvcat(sectstrmat,line);
          end
        end
      end
      [line, tmpsect] = strtok(tmpsect,eol);
    end

  end
end
