function [peakanalysis,centre, cencriteria, peak,xapp,yapp,shutterzones,xoap,errapp] = sp_pars(x,y,err,cmd,selected)
% [peakanalysis,centre, cencriteria, peak,xapp,yapp,shutterzones,xoap,errapp] = sp_pars(x,y,err,cmd,peakanalysis_in)

% uses : sp_gappf
% store : shutterzones, centre, peak, cencriteria
% show  : sp_type

if ~exist('cmd') | isempty(cmd)
	cmd = 'get';
end

if nargin == 1 & isstr(x)
	cmd = x;
	x = [];
	y = [];
end

if nargin == 2 & isstr(y)
	cmd = y;
	y = x;
	x = [];
end

if nargin == 3 & isstr(err)
	cmd = err;
	err = [];
end

if ~exist('y')
	y = [];
end

if ~exist('err')
	err = [];
end

if ~exist('selected') | isempty(selected)
	selected = ones(size(y));
	peakanalysis = [];
else
	peakanalysis = selected;
end

if ~exist('x') | isempty(x)
	if length(y) > 0
		x = 1:length(y);
	else
		x = [];
	end
end

if ~isempty(y) & (~exist('err') | isempty(err))
	err = sqrt(y-min(y)+min(y(find(y>0))));
end

if length(err) ~= length(y)
	err = ones(size(y));
end

x=x(:); y = y(:); err = err(:); shutterindex = [];
shutterzones = []; peak = []; centre = 512; cencriteria = 0; shutter = [];
xapp = []; yapp = []; xoap = []; errapp = [];
ly = length(y);

if findstr(cmd,'force')
	delete(findobj('tag','hsp_pars'));
end

if isempty(findobj('tag','hsp_pars'))

	Fposl=200; Fposb=200; Fwidth=300; Fheight=140;
	hsp_pars = figure('Position',[Fposl Fposb Fwidth Fheight],...
		'Color',get(0,'DefaultUicontrolBackgroundColor'), ...
		'MenuBar','none',...
		'Name','SSym : Parameters',...
		'NumberTitle','off',...
		'Resize','off',...
                'Tag','hsp_pars', ...
		'Visible','off');

	uicontrol(hsp_pars,...		% centre
		'String','Centre',...
		'FontName','Times','FontSize',12,...
		'Position',[5 54 50 22],...
		'Style','text',...
		'ForegroundColor',[1 1 1],...
		'HorizontalAlignment','left');

	
	uicontrol(hsp_pars,...
		'Tag','hsp_pcent',...
		'Position',[95 54 50 22],...
		'Style','edit',...
		'String','512',...
		'ForegroundColor',[0 0 0],...
		'BackgroundColor',[1 1 1],...
		'HorizontalAlignment','right');


	uicontrol(hsp_pars,...		% CenCriteria
		'String','Sym. Criteria',...
		'FontName','Times','FontSize',12,...
		'Position',[5 90 90 22],...
		'Style','text',...
		'ForegroundColor',[1 1 1],...
		'HorizontalAlignment','left');

	
	uicontrol(hsp_pars,...
		'Tag','hsp_pcrit',...
		'Position',[95 90 50 22],...
		'Style','edit',...
		'String','0',...
		'ForegroundColor',[0 0 0],...
		'BackgroundColor',[1 1 1],...
		'HorizontalAlignment','right');

	uicontrol(hsp_pars,...		% ShutterZones
		'String','Apparatus Zones',...
		'FontName','Times','FontSize',12,...
		'Position',[160 105 150 22],...
		'Style','text',...
		'ForegroundColor',[1 1 1],...
		'HorizontalAlignment','left');
	
	uicontrol(hsp_pars,...
		'Tag','hsp_pshut',...
		'Position',[160 10 130 100],...
		'Style','edit',...
		'Min',0,'Max',2,...
		'String','',...
		'ForegroundColor',[0 0 0],...
		'BackgroundColor',[1 1 1],...
		'HorizontalAlignment','right');

	uicontrol(hsp_pars,...
		'Tag','hsp_pshutind',...
		'userdata',[],...
		'Visible','off');

	uicontrol(hsp_pars,...
		'Tag','hsp_ppeak',...
		'Position',[10 150 500 100],...
		'Min',0,'Max',2,...
		'Style','text',...
		'String','',...
		'ForegroundColor',[0 0 0],...
		'BackgroundColor',[1 1 1],...
		'Visible','off');

% buttons : OK(hide), cancel

	uicontrol(hsp_pars,...
		'Style','PushButton',...
		'String','Ok',...
		'Callback','sp_pars(''hide'');',...
		'Position',[32 20 50 20]);
   	uicontrol(hsp_pars,...
		'Style','PushButton',...
		'String','Cancel',...
		'CallBack','delete(gcf)', ...
		'Position',[92 20 50 20]);

else
	hsp_pars = findobj('tag','hsp_pars');
end

% now get data from uiobjects ----------------------------------------

h=findobj(hsp_pars,'tag','hsp_ppeak');
if ~isempty(h)
	peak = str2num(get(h,'string'));
end

h=findobj(hsp_pars,'tag','hsp_pshutind');
if ~isempty(h)
	shutterindex = get(h,'userdata');
end

h=findobj(hsp_pars,'tag','hsp_pshut');
if ~isempty(h)
	shutterzones = str2num(str2mat('[ ',get(h,'string'),' ]'));
end

h=findobj(hsp_pars,'tag','hsp_pcent');
if ~isempty(h)
	centre = str2num(get(h,'string'));
end

h=findobj(hsp_pars,'tag','hsp_pcrit');
if ~isempty(h)
	cencriteria = str2num(get(h,'string'));
end

% now execute commands ----------------------------------------

if findstr(cmd,'hide')
	set(hsp_pars,'Visible','off');
	ssym;	% update plot
end

if findstr(cmd,'show')
	set(hsp_pars,'Visible','on');
	figure(hsp_pars);
end

if (~isempty(findstr(cmd,'compute')) | isempty(peak) | (isempty(shutterindex) & isempty(selected))| isempty(shutterzones)) & ~isempty(y)

	[shutterindex, peak, peakanalysis, levely, sp_type, pmin,pmax] = sp_gappf(x,y,[],[],[],(selected > 0));
	fprintf(1,'Method used was : %s\n',sp_type);
	[centre, cencriteria] = sp_cent(x,y,peak);
end

if ~isempty(shutterindex) & ~isempty(x) & (~isempty(findstr(cmd,'compute')) | isempty(shutterzones) | isempty(peakanalysis))
	shutterzones = x(shutterindex);
	[n,p] = size(shutterindex);
	shutterzones = reshape(shutterzones,n,p);
end

if ~isempty(y) & ~isempty(shutterzones) & isempty(peakanalysis)

	shutter = 0*y;
	for i=1:size(shutterzones)
		b1 = shutterzones(i,1);
		b1 = find(x >= b1);
		if isempty(b1)
			b1 = length(y);
		else
			b1 = b1(1);
		end
		b2 = shutterzones(i,2);
		b2 = find(x <= b2);
		if isempty(b2)
			b2 = 1;
		else
			b2 = b2(length(b2));
		end
		shutter(b1:b2) = 1;
	end

% separate data into background, apparatus and peak parts.

	np = size(peak,1); % number of significant peaks
	peakanalysis = 3*ones(size(x));
	for i = 1:np
		peakanalysis(max(1,peak(i,1)-2*peak(i,2)):min(ly,peak(i,1)+2*peak(i,3))) = 2;
	end
	peakanalysis(find(shutter)) = 1;
	i = find(selected <= 0);
	peakanalysis(i) = -abs(peakanalysis(i));
end
if ~isempty(peak) & ~isempty(peakanalysis)

	peakanalysis = peakanalysis(:);
	% get yapp data from peak matrix

	h = findobj('tag','hsp_ogapp');
	oglobapp = get(h,'Value') - 1;

	i = find(peakanalysis(peak(:,1)) == 1);	% indexes of apparatus peaks
	if ~isempty(i)

	peakapp = peak(i,:);

	y = y.*(peakanalysis == 1);
	int = y(peakapp(:,1))-min(y(peakapp(:,7)),y(peakapp(:,8)));

	[dummy,yapp] = max(int); % biggest thinner peak in apparatus zones
	xoap = x(peakapp(yapp,1));
	pos = peakapp(yapp,1);	% working on indexes
	width = 3*max(peakapp(yapp,2),peakapp(yapp,3));
	% yapp = y( peakapp(yapp,7):peakapp(yapp,8) );
	if pos-width < 1
		width = pos -1;
	end
	if pos + width > ly
		width = (ly-pos);
	end
	id = (pos-width):(pos+width);
	yapp = y( id ).*(peakanalysis( id ) == 1); % odd number of point for exact centering
	xapp = x( id );
	errapp = err( id );

	[int,pos] = max(yapp);	% now centering
	width = min(pos -1, length(yapp) - pos);
	id = (pos-width):(pos+width);
	yapp = yapp(id);
	errapp = errapp(id);
	xapp=xapp(id);
	yapp = yapp-min(yapp);
	yapp (1) = 0;
	yapp(end) = 0;

	yapp = yapp/sum(yapp);

	end % if ~empty(i)
end

% now set data to uiobjects ----------------------------------------

h=findobj(hsp_pars,'tag','hsp_ppeak');
if ~isempty(h)
	set(h,'string',num2str(peak));
	if isempty(str2num(get(h,'string')))
		set(h,'string',num2str(peak,'%.8g '));
	end
end

h=findobj(hsp_pars,'tag','hsp_pshut');
if ~isempty(h)
	set(h,'string',num2str(shutterzones,'%.8g '));
end

h=findobj(hsp_pars,'tag','hsp_pshutind');
if ~isempty(h)
	set(h,'userdata',shutterindex);
end
h=findobj(hsp_pars,'tag','hsp_pcent');
if ~isempty(h)
	set(h,'string',num2str(centre));
end

h=findobj(hsp_pars,'tag','hsp_pcrit');
if ~isempty(h)
	set(h,'string',num2str(cencriteria));
end

peakanalysis = peakanalysis(:);
