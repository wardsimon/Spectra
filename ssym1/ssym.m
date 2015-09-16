function [x,y,err,selected,peakanalysis,xapp,yapp] = ssym(x,y,err,selected,cmd)
% SSym     : Spectrum symetrization, with MFit
% [x,y,err,selected,peakanalysis,xapp,yapp] = ssym(x,y,err,selected,cmd) or ssym(cmd)
% Note : MFit is needed (mf_flsqr, mf_pwin, mf_fndpks...)
% cmd :
%		force
%		exit
%		symx
%		symy
%		tomfit
%		fromfit
% 		togglelinlog
%		togglegrid
%		auto{x|y|c} on|off
%		about

ssym_ver = '0.1';

markersize = 6;
marker = '.';
bgcolor = 'cyan';
sigcolor = 'magenta';
appcolor = 'green';
uselcolor = 'blue';

if ~exist('cmd') | isempty(cmd)
	cmd = 'update';
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

if nargin == 4 & isstr(selected)
	cmd = selected;
	selected = [];
end

if ~exist('y')
	y = [];
end

if ~exist('err')
	err = [];
end

if ~exist('x')
	x = [];
end

if ~exist('selected')
	selected = [];
end

x=x(:); y = y(:); err = err(:); selected = selected(:);
xapp = []; yapp = []; peakanalysis = []; xoap = []; errapp = [];
order=[]; optionx=[]; appfcn=[]; specfcn=[]; optiony=[];

name = '';

if findstr(cmd,'force')
	delete(findobj('tag','hsp_ctrl'));
end

if isempty(findobj('tag','hsp_ctrl'))
	if ~exist('mfit')
		disp('Sorry this program needs MFit package to run')
		return
	end
	disp([ 'Welcome to SSym ' ssym_ver ]);

	Fposl=10; Fposb=200; Fwidth=460; Fheight=400;
	hsp_ctrl = figure('Position',[Fposl Fposb Fwidth Fheight],...
		'Color','black', ...
		'MenuBar','figure',...
		'Name','SSym : Control',...
		'NumberTitle','off',...
		'Resize','on',...
                'Tag','hsp_ctrl', ...
		'HandleVisibility','on',...
		'DefaultUicontrolHorizontalAlignment','center',...
		'DeleteFcn','ssym(''exit'');',...
		'WindowButtonDownFcn','sp_rbbox(''go'');',...
		'Visible','on');
	uicontrol(hsp_ctrl,'Tag','sp_zoom_list',...
		'visible','off');      % Zoom list
	uicontrol(hsp_ctrl,'Tag','sp_rbbox_action',...
		'visible','off','userdata','zoomin');   % Action of rbbox

	uicontrol(hsp_ctrl,'Tag','sp_xoap',...
		'visible','off','userdata',512);   % Centre of App in (x,y)

	hsp_plot = axes('position',[ 0.1 0.55 0.8 0.4 ],'tag','hsp_plot');
	hsp_app  = axes('position',[ 0.1 0.05 0.4 0.4 ],'tag','hsp_app');

	hsp_msym = uimenu(hsp_ctrl,'Label','&Symetrization','tag','hsp_msym');
	uimenu(hsp_msym,'Label','X axis',...
		'accelerator','x',...
		'CallBack','ssym(''symx'');');
	uimenu(hsp_msym,'Label','Y axis',...
		'accelerator','y',...
		'CallBack','ssym(''symy'');');
	uimenu(hsp_msym,'Label','Change X axis...',...
		'CallBack','mf_chgx(''ssym+show'');');
	uimenu(hsp_msym,'Label','Options...',...
		'Separator','on',...
		'CallBack','sp_opt(''show'');');
	uimenu(hsp_msym,'Label','Parameters...',...
		'CallBack','sp_pars(''show'');');
	uimenu(hsp_msym,'Label','Apparatus Fcn...',...
		'CallBack','spapp(''show'');');


	hsp_mdat = uimenu(hsp_ctrl,'Label','&Data','tag','hsp_mdat');
	uimenu(hsp_mdat,'Label','Get From Mfit',...
			'Accelerator','g',...
			'CallBack','ssym(''fromfit'');');
	uimenu(hsp_mdat,'Label','Send To Mfit',...
			'Accelerator','s',...
			'CallBack','ssym(''tomfit'');');
	uimenu(hsp_mdat,'Label','Get From Matlab','CallBack','ssym(''frombase'');');
	uimenu(hsp_mdat,'Label','Transform...','CallBack','ssym(''fromssym'');');
	uimenu(hsp_mdat,...
			'Label','Select Data',...
			'Separator','on',...
			'Tag','hsp_mvsel',...
			'Checked','off',...
			'Callback','sp_rbbox(''select'')');
	uimenu(hsp_mdat,...
			'Label','Select Peak',...
			'Tag','hsp_mvsep',...
			'Checked','off',...
			'Accelerator','k',...
			'ForegroundColor',sigcolor,...
			'Callback','sp_rbbox(''selectpk'')');
	uimenu(hsp_mdat,...
			'Label','Select Background',...
			'Tag','hsp_mvseb',...
			'Checked','off',...
			'ForegroundColor',bgcolor,...
			'Accelerator','b',...
			'Callback','sp_rbbox(''selectbg'')');
	uimenu(hsp_mdat,...
			'Label','Select Apparatus',...
			'accelerator','a',...
			'Tag','hsp_mvsea',...
			'Checked','off',...
			'ForegroundColor',appcolor,...
			'Callback','sp_rbbox(''selectapp'')');
	uimenu(hsp_mdat,...
			'Label','Deselect Data',...
			'Tag','hsp_mvusel',...
			'Checked','off',...
			'ForegroundColor',uselcolor,...
			'Callback','sp_rbbox(''deselect'')');
	uimenu(hsp_mdat,...
			'Label','Select All',...
			'Callback','sp_rbbox(''selectall'')');
	uimenu(hsp_mdat,...
			'Label','Deselect All',...
			'Callback','sp_rbbox(''deselectall'')');
	uimenu(hsp_mdat,'Label','About...',...
		'Separator','on',...
		'CallBack','ssym(''about'');');
	uimenu(hsp_mdat,'Label','Exit',...
		'accelerator','q',...
		'CallBack','ssym(''exit'');');

	hsp_mvue = uimenu(hsp_ctrl,'Label','&View','tag','hsp_mvue');
	uimenu(hsp_mvue,...
		'Tag','hsp_mvllog',...
		'Label','Lin/Log Y axis',...
		'accelerator','l',...
		'Checked','off',...
		'CallBack','ssym(''togglelinlog'');');
	uimenu(hsp_mvue,...
		'Tag','hsp_mvgrid',...
		'Label','Grid On/Off',...
		'Checked','off',...
		'CallBack','ssym(''togglegrid'');');
	uimenu(hsp_mvue,...
		'Tag','hsp_mvzoom',...
		'Label','Zoom In',...
		'Checked','on',...
		'Accelerator','z',...
		'CallBack','sp_rbbox(''zoomin'');');
	uimenu(hsp_mvue,...
		'Label','Zoom Out',...
		'CallBack','sp_rbbox(''zoomout'');');
	uimenu(hsp_mvue,...
		'Label','Reset Zoom',...
		'Separator','on',...
		'Callback','sp_rbbox(''zoomreset'')');
	uimenu(hsp_mvue,...
		'Label','Replot',...
		'accelerator','r',...
		'Callback','ssym;');
	uimenu(hsp_mvue,...
		'Label','Rescan',...
		'Callback','ssym(''init'');');

% now some 'quick' buttons and auto-mode settings

	uicontrol(hsp_ctrl,...
		'Units','normalized',...
		'Style','Frame',...
		'BackgroundColor',[0 0 0],...
		'Position',[ 0.51 0.05 0.48 0.35 ]);

	uicontrol(hsp_ctrl,...
		'Units','normalized',...
		'Style','PushButton',...
		'HorizontalAlignment','center',...
		'String','Sym X',...
		'Callback','ssym(''symx'');',...
		'Position',[ 0.7 0.21 0.2 0.06 ]);
	uicontrol(hsp_ctrl,...
		'Units','normalized',...
		'Style','PushButton',...
		'HorizontalAlignment','center',...
		'String','Sym Y',...
		'Callback','ssym(''symy'');',...
		'Position',[ 0.7 0.14 0.2 0.06]);
	uicontrol(hsp_ctrl,...
		'Units','normalized',...
		'Style','PushButton',...
		'HorizontalAlignment','center',...
		'String','Chg X axis...',...
		'Callback','mf_chgx(''ssym+show'');',...
		'Position',[ 0.7 0.07 0.2 0.06]);
	uicontrol(hsp_ctrl,...
		'Units','normalized',...
		'Style','checkbox',....
		'Tag','hsp_bautoy',...
		'String','AutoY',...
		'Position',[ 0.57 0.15 0.12 0.04]);
	uicontrol(hsp_ctrl,...
		'Units','normalized',...
		'Style','checkbox',....
		'Tag','hsp_bautox',...
		'String','AutoX',...
		'Position',[ 0.57 0.22 0.12 0.04]);

	uicontrol(hsp_ctrl,...
		'Units','normalized',...
		'Style','checkbox',....
		'Tag','hsp_bautoc',...
		'String','Auto',...
		'Position',[ 0.57 0.08 0.12 0.04]);

	uicontrol(hsp_ctrl,...
		'Units','normalized',...
		'Style','PushButton',...
		'HorizontalAlignment','center',...
		'String','Params...',...
		'Callback','sp_pars(''show'');',...
		'Position',[ 0.54 0.30 0.14 0.06]);

	uicontrol(hsp_ctrl,...
		'Units','normalized',...
		'Style','PushButton',...
		'String','Options...',...
		'Callback','sp_opt(''show'');',...
		'Position',[ 0.69 0.30 0.14 0.06]);

	uicontrol(hsp_ctrl,...
		'Units','normalized',...
		'Style','PushButton',...
		'String','Yapp...',...
		'Callback','spapp(''show'');',...
		'Position',[ 0.84 0.30 0.14 0.06]);

	set(hsp_ctrl,'userdata',[ hsp_plot hsp_app ]);
else
	hsp_ctrl = findobj('tag','hsp_ctrl');
	userdata = get(hsp_ctrl,'userdata');
	if ~isempty(userdata)
		hsp_plot = userdata(1);
		hsp_app  = userdata(2);
	end
end

% now import data --------------------------------------------------

if findstr(cmd,'fromfit')
	[x,y,err,selected]=fromfit;
	h = findobj('Tag','mf_DataFile');
	if ~isempty(h)
		name = get(h,'string');
	end
end

if findstr(cmd,'frombase')
	[x, y, err]=frombase('base');
	selected = [];
	name = 'User';
end

parcmd = 'get';

if ~isempty(y) | ~isempty(x)
	parcmd = 'compute';	% new data
end

if isempty(y) | isempty(x) | isempty(err) | isempty(selected) | isempty(peakanalysis)
	userdata = get(hsp_plot,'userdata');
	if ~isempty(userdata)
		if isempty(x), x = userdata(:,1); end
		if isempty(y), y = userdata(:,2); end
		if isempty(err), err = userdata(:,3); end
		if isempty(selected), selected = userdata(:,4); end
		if isempty(peakanalysis), peakanalysis = userdata(:,5); end
	end
	userdata = get(hsp_app,'userdata');
	if ~isempty(userdata)
		yapp = userdata(:,2);
		xapp = userdata(:,1); % should be centered
	end
end

if findstr(cmd,'fromssym')
	[x, y, err]=frombase(x,y,err,'SSym :');
	selected = [];
	name = 'SSym';
end

if length(x) ~= length(y) & length(y) > 0
	x = 1:length(y);
end

if ~isempty(y) & isempty(err)
	err = sqrt(y);
	i = find(err == 0);
	err(i) = sqrt(min(y(find(y>0))));
end

if length(err) ~= length(y)
	err = ones(size(y));
end

i = find(y <= 0);
if ~isempty(i)
	y = y - min(y);
	i = find(y <= 0);
	y(i) = min(y(find(y>0))/10);
end

if length(selected) ~= length(y)
	selected = ones(size(y));
end

if length(peakanalysis) ~= length(y)
	peakanalysis = [];
end

% init commands ----------------------------------------

if findstr(cmd,'autox')
	if findstr(cmd,'on')
		set(findobj('Tag','hsp_bautox'),'Value',1);
	else
		set(findobj('Tag','hsp_bautox'),'Value',0);
	end
end

if findstr(cmd,'autoy')
	if findstr(cmd,'on')
		set(findobj('Tag','hsp_bautoy'),'Value',1);
	else
		set(findobj('Tag','hsp_bautoy'),'Value',0);
	end
end

if findstr(cmd,'autoc')
	if findstr(cmd,'on')
		set(findobj('Tag','hsp_bautoc'),'Value',1);
	else
		set(findobj('Tag','hsp_bautoc'),'Value',0);
	end
end


if ~isempty(findstr(cmd,'init')) & isempty(findstr(cmd,'exit'))	% force computation, normally done on loading
	parcmd = 'compute';
	peakanalysis = [];
end

[peakanalysis, centre, cencriteria, peak, xapp, yapp, shutterzones, xoap, errapp] = sp_pars(x,y,err,parcmd,peakanalysis);

[order, optionx, appfcn, specfcn, optiony] = sp_opt('get');

% 'compute' commands ----------------------------------------

i = find(~selected);
peakanalysis(i) = -abs(peakanalysis(i));	% do not take into account

autox = get(findobj('Tag','hsp_bautox'),'Value');
if ~isempty(findstr(cmd,'symx')) | autox
	[newx, centre, peak, peakanalysis, rescaletable] = sp_symx(x, y, err, order, optionx, peakanalysis, centre, cencriteria, peak, appfcn, specfcn);
	if isempty(rescaletable)
		disp('Nothing was done on X')
	else
	% transfert to sp_pars : peak, centre

	hsp_pars = findobj('tag','hsp_pars');
	h=findobj(hsp_pars,'tag','hsp_ppeak');
	if ~isempty(h)
		set(h,'string',num2str(peak));
	else
		disp('Can''t set Peak matrix')
	end

	h=findobj(hsp_pars,'tag','hsp_pcent');
	if ~isempty(h)
		set(h,'string',num2str(centre));
	else
		disp('Can''t set Centre')
	end
	[x, sorti]  = sort(newx);
	y = y(sorti);
	disp('X sym done.')
	[peakanalysis, centre, cencriteria, peak, xapp, yapp,shutterzones,xoap,errapp] = sp_pars(x,y,err,'update',peakanalysis);
	parcmd = 'compute';
	end
end

autoy = get(findobj('Tag','hsp_bautoy'),'Value');
if ~isempty(findstr(cmd,'symy')) | autoy
	if isempty(findstr(cmd,'symx')) | ~isempty(rescaletable)
		[newx,newy,newerr,peakanalysis] = sp_symy(x,y,err,centre,peak,peakanalysis,[optiony order]);
		x = newx; y = newy;
		disp('Y sym done.')
		[peakanalysis, centre, cencriteria, peak, xapp, yapp,shutterzones,xoap,errapp] = sp_pars(x,y,err,'update',peakanalysis);
		parcmd = 'compute';
	else
		disp('Nothing was done on Y');
	end
end

i = find(~selected);
peakanalysis(i) = -abs(peakanalysis(i));	% do not take into account

% update all graphs and figures ----------------------------------------
x=x(:); y = y(:); err = err(:);
peakanalysis = peakanalysis(:); selected = selected(:);
xapp = xapp(:); yapp = yapp(:); errapp = errapp(:);

if (~isempty(x) & ~isempty(y))	% put data into ssym plot 
	axes(hsp_plot);
	curax = [get(hsp_plot,'XLim') get(hsp_plot,'YLim')];
	linlog = get(hsp_plot,'Yscale');
	i = find(x>curax(1) & x<curax(2) & y>curax(3) & y<curax(4));
	if ~isempty(i)
		markersize2 = min(20,markersize*length(x)/length(i));
	else
		markersize2 = markersize;
	end
	hold off
	cla
	i = find(peakanalysis==3 & y > 0);	% background
	plot(x(i),y(i),'Marker',marker,...
		'LineStyle','none',...
		'Tag','hsp_plsel',...
		'MarkerSize',markersize2,'Color',bgcolor);
	hold on
	i = find(peakanalysis==1 & y > 0);	% apparatus
	plot(x(i),y(i),'Marker',marker,...
		'LineStyle','none',...
		'Tag','hsp_plsel',...
		'MarkerSize',markersize2,'Color',appcolor);
	i = find(peakanalysis==2 & y > 0);	% signal
	plot(x(i),y(i),'Marker',marker,...
		'LineStyle','none',...
		'Tag','hsp_plsel',...
		'MarkerSize',markersize2,'Color',sigcolor);
	i = find(peakanalysis<=0 & y > 0);	% unselected
	plot(x(i),y(i),'Marker',marker,...
		'LineStyle','none',...
		'MarkerSize',ceil(markersize/2),'Color',uselcolor);
	set(hsp_plot,'Xlim',curax(1:2));
	set(hsp_plot,'Ylim',curax(3:4),'Yscale',linlog);
	hold off
	ylabel('signal')
end

if ~isempty(x) & ~isempty(y) & ~isempty(err) & ~isempty(peakanalysis)
	set(hsp_plot,'userdata',[ x y err selected peakanalysis ]);
end

if (~isempty(xapp) & ~isempty(yapp))	% put data into ssym plot app
	axes(hsp_app);
	axis auto
	h = findobj(hsp_app,'Tag','hsp_fitapp');
	if ~isempty(h)
		hold on
	end
	h = findobj(hsp_app,'Tag','hsp_plapp');
	if ~isempty(h)
		delete(h)
	end
	h = plot(xapp,yapp);
	hold off
	set(h,'Color',appcolor,'tag','hsp_plapp');
	tmp = axis;
	tmp(1) = min(xapp); tmp(2) = max(xapp);
	axis(tmp);
	t = 'Apparatus function (spapp)';
	if ~isempty(xoap)
		t = [ t ' at x=' num2str(xoap,3) ];
	end
	title(t);
end

if ~isempty(xapp) & ~isempty(yapp)
	set(hsp_app,'userdata',[ xapp yapp errapp ]);
end

if ~isempty(xoap)
	set(findobj('Tag','sp_xoap'),'userdata',xoap);
end

% now doing minor commands ----------------------------------------

if findstr(cmd,'togglelinlog')
	if strcmp(get(findobj('Tag','hsp_mvllog'),'Checked'),'off')
		set(findobj('Tag','hsp_mvllog'),'Checked','on');
		set(hsp_plot,'yscale','log');
		set(hsp_app,'yscale','log');
	else
		set(findobj('Tag','hsp_mvllog'),'Checked','off');
		set(hsp_plot,'yscale','linear');
		set(hsp_app,'yscale','linear');
	end
	axes(hsp_app);
	fcns = get(findobj('Tag','hsp_mvgrid'),'Checked');
	axis auto
	tmp = axis;
	if ~isempty(tmp) & ~isempty(xapp)
		tmp(1) = min(xapp); tmp(2) = max(xapp);
		axis(tmp);
	end
	set(hsp_plot,'Xgrid',fcns,'Ygrid',fcns);
	set(hsp_app,'Xgrid',fcns,'Ygrid',fcns);
end

if findstr(cmd,'togglegrid')
	fcns = get(findobj('Tag','hsp_mvgrid'),'Checked');
	if findstr(fcns,'off')
		fcns = 'on';
	else
		fcns = 'off';
	end
	set(hsp_plot,'Xgrid',fcns,'Ygrid',fcns);
	set(hsp_app,'Xgrid',fcns,'Ygrid',fcns);
	set(findobj('Tag','hsp_mvgrid'),'Checked',fcns);
	% drawnow;
end

if findstr(cmd,'about')
	helpstr = {'SSym : This program enables to',...
	'1- get data from MFit or Matlab',...
	'2- view/zoom/select data',...
	'3- search for symetry',...
	'4- force symetry on X and Y axis',...
	'5- change X axis',...
	'6- configure many features at users choice',...
	' ',...
	'Authors : EF (1998)'};
	helpdlg(helpstr,'SSym About');
end

% ending commands ----------------------------------------

if findstr(cmd,'tomfit')
	tomfit(x,y,err,(selected & (peakanalysis ~= 1)));
	parcmd = 'compute';
end

autoc = get(findobj('Tag','hsp_bautoc'),'Value');
if autoc & strcmp(parcmd,'compute') & isempty(findstr(cmd,'init'))
	mf_chgx('ssym+newx');
	if nargout <= 1
		x = [ cmd ' done' ];
	end
	return;
end

if findstr(cmd,'exit')
	%===== Close the windows ============================================
	hmat = get(0,'children');
	for i = 1:size(hmat)
		h = hmat(i);
		tag = eval('get(h,''tag'')','[]');
		if ~isempty(strmatch('sp',tag)) | ~isempty(strmatch('hsp',tag))
			delete(h);
		end
	end
end

if strcmp(parcmd,'compute')
	sp_rbbox('zoomreset');
	set(findobj('Tag','hsp_mvllog'),'Checked','off');
	set(findobj('Tag','hsp_mvgrid'),'Checked','off');
end

if nargout <= 1
	x = [ cmd ' done' ];
end
