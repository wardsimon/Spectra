function makeeps(figm,position,option)
% makeeps  : prepares a figure for EPS printing
% makeeps(fig,size,option)
% size can be a scale factor or direct size in figure unit.
% option can be 'silent', 'normal, 'verbose'

if nargin < 1, figm=[]; end
if nargin < 2, position=[]; end
if nargin < 3, option=[]; end
if isempty(figm)
  figm = get(0,'CurrentFigure');
end

for fig = figm

figure(fig)

figpos = get(fig,'paperposition');   % initial size

if isempty(option)
  option = 'normal';
end

if ~isempty(findstr(option,'silent'))
  output = 0;
elseif ~isempty(findstr(option,'verbose'))
  output = 2;
else
  output = 1;
end

figpos = figpos(3:4);
if output
fprintf(1,'Initial size of figure %i is [%f x %f] %s\n',fig,figpos(1),figpos(2),get(fig,'paperunit'));
end

if isempty(position)
  position = figpos;
end

if length(position) == 1
  position = position * figpos;
end

xsc = position(1)/figpos(1);
ysc = position(2)/figpos(2);
tsc = min(xsc,ysc);
htext=findobj(fig,'type','text');

for i=1:length(htext)
  h = htext(i);
  fs = get(h,'FontSize');
  t = get(h,'string');
  if iscellstr(t), t=[ t{1} '...' ]; end
  if output > 1
  fprintf(1,'Changing text %f (%s) size %i -> %i \n',h,t,fs,ceil(fs*sqrt(tsc)));
  end
  set(h,'FontSize',ceil(fs*tsc));
end

hlin=findobj(fig,'type','line');
for i=1:length(hlin)
  h = hlin(i);
  fs = get(h,'MarkerSize');
  if output > 1
  fprintf(1,'Changing Marker %f (%s) size %i -> %i \n',h,get(h,'Marker'),fs,ceil(fs*sqrt(tsc)));
  end
  set(h,'MarkerSize',ceil(fs*tsc));
end

haxmat = findobj(fig,'type','axes');
for i = 1:length(haxmat)
  hax = haxmat(i);
  if output > 1
  fprintf(1,'Handling Axis %f\n',hax);
  end
  hlx = get(hax,'XLabel');
  hly = get(hax,'YLabel');
  hlz = get(hax,'ZLabel');
  htl = get(hax,'Title');
  for h= [ hlx hly hlz htl ]
    fs = get(h,'FontSize');
    t = get(h,'string');
    if iscellstr(t), t=[ t{1} '...' ]; end
    if output > 1
    fprintf(1,'Changing text %f (%s) size %i -> %i \n',h,t,fs,ceil(fs*sqrt(tsc)));
    end
    set(h,'FontSize',ceil(fs*sqrt(tsc)));
  end
  fs = get(hax,'FontSize');
  set(hax,'FontSize',ceil(fs*sqrt(tsc)));
end

set(fig,'PaperPositionMode','manual');
set(fig,'PaperPosition',[ 0 0 position ]);
figxy = get(fig,'position');
set(fig,'position',[ figxy(1:2) xsc*figxy(3) ysc*figxy(4) ]);
if output
fprintf(1,'Final size is [%f x %f] %s\n',position(1),position(2),get(fig,'paperunit'));
end
if output > 1
fprintf(1,'Type in ''print -f%i -loose -depsc <output>'' to print figure.\n',fig);
end
end
