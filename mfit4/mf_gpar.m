function [x,y,err,xname,yname,monitor,sectstrmat,fullsecstr] = mf_gpar(str,parnamelistx,parnamelisty,fullsecstr)
% [x,y,err,xname,yname,monitor,sectstrmat,fullsecstr ] = mf_gpar(str or file,xname,yname,sectionselect or sectionstr)
% to extract n parameters for a Par section into x and y variables
% err contains dy value if available.
% sectionselect can be a vector, 'all', 'param', or 'vars'.
% exemple :
%   [T,laser,dummy,xname,yname,monitor,str] = mf_gpar('s100pur.mft','f3:Temperature','f2:Laser','all');
%   [x,y,dy] = mf_gpar('s100pur.mft');   % works fine, and ask user....

if nargin < 2
  parnamelistx = '';
end
if nargin < 3
  parnamelisty = '';
end
if nargin < 4
  fullsecstr = [];
end
sectionselect = [];
SectionType = 'Param';
if findstr(fullsecstr,'param')
  SectionType = 'Param';
end
if findstr(fullsecstr,'vars')
  SectionType = 'Vars';
end
if ~isempty(findstr(fullsecstr,'all')) | ~ischar(fullsecstr)
  sectionselect = fullsecstr;
  fullsecstr = [];
end

x = [];
y = [];
err=[];
xname = [];
yname = [];
sectstrmat = [];
selectx = [];
selecty = [];
if size(str,1) == 1 & isempty(fullsecstr)
  fprintf(1,'Collecting %s sections in data\n', SectionType);
  if strcmp(SectionType,'Vars')
    [fullsection, sectstrmat, fullsecstr] = mft_get(str,SectionType,sectionselect,'=');
  else
    [fullsection, sectstrmat, fullsecstr] = mft_get(str,SectionType,sectionselect,'par');

  end
  % collect Parameter sections (and ask user to select sections)
  fprintf(1,'Found %i sections\n',size(fullsecstr,1));
else
  sectstrmat = str;
end


cursect = 0;
parnamelist = '';
for i=1:size(sectstrmat,1)
  line = sectstrmat(i,:);
  findsect = findstr(line,'Section');
  if ~isempty(findsect)
    cursect = i;
  end
  if cursect & strcmp(SectionType,'Param') & ~isempty(findstr(line,'par'))
    parname = sscanf(line,'%*s%*d%s%*f%*f');
    if isempty(strmatch(parname,parnamelist))
      parnamelist = strvcat(parnamelist,parname);
    end
%   pardat = sscanf(line,'%*s %d %*s %f %f %f');
  end
        if cursect & strcmp(SectionType,'Vars') & findstr(line,' = ')
                parname = sscanf(line,'%s=%*s');
                if isempty(strmatch(parname,parnamelist))
      parnamelist = strvcat(parnamelist,parname);
    end
        end
end
parnamelist = strvcat(parnamelist,sprintf('Index (1:%i)',size(fullsecstr,1)));

% now ask to select parameters to import
if isempty(parnamelistx)
  [selectx,okx] = listdlg('ListString',cellstr(parnamelist),...
  'PromptString','Select X parameters to import',...
  'Name','Set X Parameters','ListSize',[300 160]);
  if ~okx
    return
  end
  parnamelistx = parnamelist(selectx,:);
end


if isempty(parnamelisty)
  [selecty,oky] = listdlg('ListString',cellstr(parnamelist),...
  'PromptString','Select Y parameters to import',...
  'Name','Set Y Parameters','ListSize',[300 160]);
  if ~oky
    return
  end
  parnamelisty = parnamelist(selecty,:);
end

pardatx = [];
pardaty = [];
cursec = 0;
secindex = 0;

for i=1:size(sectstrmat,1)
  line = sectstrmat(i,:);

  if cursect & ~isempty(findstr(line,'par')) & strcmp(SectionType,'Param')
    parname = sscanf(line,'%*s%*d%s%*f%*f');
    dat = sscanf(line,'%*s %d %*s %f %f %f');
    if length(dat) == 3
      dat = [ dat(:) ; 0];
    end

    if ~isempty(strmatch(parname,parnamelistx))
      pardatx = [ pardatx dat(:) ];
    end
    if ~isempty(strmatch(parname,parnamelisty))
      pardaty = [ pardaty dat(:) ];
    end
    if ~isempty(findstr(parnamelistx,sprintf('Index (1:%i)',size(fullsecstr,1))))
      pardatx = [ 0 secindex 0 0 ]';
    end
  end
  oky = findstr(line,'=');
  if cursect & ~isempty(oky) & strcmp(SectionType,'Vars')
    parname = sscanf(line,'%s=%*s');
    dat = line((oky(1)+1):end);

    if ~isempty(strmatch(parname,parnamelistx))
      pardatx = [ pardatx dat ];
    end
    if ~isempty(strmatch(parname,parnamelisty))
      pardaty = [ pardaty dat ];
    end
  end

  findsect = findstr(line,'Section'); % new section begins or EOF ?
  if ~isempty(findsect) | (i == size(sectstrmat,1))
    if i < size(sectstrmat,1)
      cursect = i;
      secindex = secindex +1;
    end
    if ischar(pardatx)
      x = [ x ; pardatx ];
      y = [ y ; pardaty ];
    else
      if size(pardatx,2) == 1 & size(pardaty,2) > 1
        pardatx = kron(pardatx,ones(1,size(pardaty,2)));
      elseif size(pardaty,2) == 1 & size(pardatx,2) > 1
        pardaty = kron(pardaty,ones(1,size(pardatx,2)));
      end
      if length(pardatx) == length(pardaty)
        x = [ x pardatx ];
        y = [ y pardaty ];
      end
    end
    pardatx = [];
    pardaty = [];
  end
end
if ~isempty(x) & ~ischar(x) & (size(x,1) == size(y,1)) & (size(x,2) == size(y,2))
  [dummy, sorti] = sort(x(1,:));
  x = x(2,sorti);
  err = y(3,sorti);
  y = y(2,sorti);
end
if ~isempty(selectx)
  xname = parnamelist(selectx(1),:);
else
  xname = parnamelistx;
end
if ~isempty(selecty)
  yname = parnamelist(selecty(1),:);
else
  yname = parnamelisty;
end

xname = sprintf('%s (%i sec.)',xname,size(fullsecstr,1));
monitor = ones(size(y));
