function [yapp, name, pnames, pin] = spapp(x,pin,flag)
% [y, name, pnames, pin] = spapp(x,p,flag)
% This is the apparatus function determined by SSym
% Normal function with x,p
% according to Exp/Model flag, use experimental
% data (no pars) or fit by a function (appfcn).
% other direct commands as usual : hide, show, get, fit(and guess from peak)

y = []; 
if ~exist('x')
	x=[];
end
yapp = 0*x;
name = 'SSym';
pnames = [];

% check for SSym ----------------------------------------

hsp_plot = []; hsp_app = [];
hsp_ctrl = findobj('tag','hsp_ctrl');
hsp_pars = findobj('tag','hsp_pars');
userdata = get(hsp_ctrl,'userdata');
if ~isempty(userdata)
	hsp_plot = userdata(1);
	hsp_app  = userdata(2);
else
	disp('SSym is not running.')
	return
end

% check parameters ----------------------------------------

cmd = ''; peak = [];

if nargin == 1 & isstr(x)
	cmd = x;
	x = [];
end

if nargin == 3
	if length(flag) == 1
		switch (flag)
			case 0, cmd = [ cmd '+get'];
			case 1, cmd = [ cmd '+identify' ];
			case 2, cmd = [ cmd '+guess' ];
			case 3, cmd = [ cmd '+guess+fit' ];
		end
	elseif length(flag) > 1
		cmd = [ cmd '+fit' ];
	end	
end

if nargin == 2
	cmd = '';	% retreive apparatus function : get mode
end

if exist('x')
	xin = x;
else
	xin = [];
end

if ~exist('pin')
	pin = [];
end

% build window if necessary

if findstr(cmd,'force')
	delete(findobj('tag','hsp_wapp'));
end

if isempty(findobj('tag','hsp_wapp'))

	Fposl=200; Fposb=350; Fwidth=300; Fheight=60;
	hsp_wapp = figure('Position',[Fposl Fposb Fwidth Fheight],...
		'Color',get(0,'DefaultUicontrolBackgroundColor'), ...
		'MenuBar','none',...
		'Name','SSym : App Fcn',...
		'NumberTitle','off',...
		'Resize','off',...
		'userdata',[],...
		'Tag','hsp_wapp', ...
		'Visible','off');

	uicontrol(hsp_wapp,...		% Uses
		'String','Uses',...
		'FontName','Times','FontSize',12,...
		'Position',[ 4 33 50 22],...
		'Style','text',...
		'ForegroundColor',[1 1 1],...
		'HorizontalAlignment','left');

	uicontrol(hsp_wapp,...		% Fit or Data
		'Style','popup',....
		'Tag','hsp_adatfit',...
		'String','Data|Fit',...
		'Value',1,...
		'Position',[ 50 35 50 20 ]);

	uicontrol(hsp_wapp,...		% With
		'String','with :',...
		'FontName','Times','FontSize',12,...
		'Position',[ 106 33 50 22],...
		'Style','text',...
		'ForegroundColor',[1 1 1],...
		'HorizontalAlignment','left');

	uicontrol(hsp_wapp,...		% app fcn
		'Style','edit',....
		'Tag','hsp_afapp',...
		'String','gauss',...
		'BackgroundColor',[1 1 1],...
		'ForegroundColor',[0 0 0],...
		'Position',[ 150 35 70 20 ]);

	uicontrol(hsp_wapp,...		% auto fit mode
		'Style','popup',....
		'Tag','hsp_afauto',...
		'String','Auto|Manual',...
		'Value',1,...
		'Position',[ 220 35 60 20 ]);


% buttons : OK(hide), cancel, Fit

	uicontrol(hsp_wapp,...
		'Style','PushButton',...
		'String','Ok',...
		'Callback','spapp(''hide'');',...
		'Position',[32 4 50 20]);
   	uicontrol(hsp_wapp,...
		'Style','PushButton',...
		'String','Cancel',...
		'CallBack','delete(gcf)', ...
		'Position',[92 4 50 20]);

   	uicontrol(hsp_wapp,...
		'Style','PushButton',...
		'String','Fit',...
		'CallBack','spapp(''fit'');', ...
		'Position',[152 4 50 20]);

   	uicontrol(hsp_wapp,...
		'Style','PushButton',...
		'String','Guess',...
		'CallBack','spapp(''guess'');', ...
		'Position',[ 212 4 50 20]);

else
	hsp_wapp = findobj('Tag','hsp_wapp');
end

% get useful data ----------------------------------------

userdata = get(hsp_app,'userdata');
if ~isempty(userdata)
	xapp = userdata(:,1);
	yapp = userdata(:,2);
	errapp = userdata(:,3);
else	
	disp('No SSym App data')
	return
end

h=findobj(hsp_wapp,'tag','hsp_afapp');
if ~isempty(h)
	appfcn = get(h,'string');
else
	disp('Can''t get AppFcn name')
	return
end

h=findobj(hsp_wapp,'tag','hsp_adatfit');
if ~isempty(h)
	datfit = get(h,'Value');	% 1=data, 2=fit
else
	disp('Can''t get DatFit mode')
	return
end

userdata = get(hsp_plot,'userdata');
if ~isempty(userdata)
	x = userdata(:,1);
	y = userdata(:,2);
	err = userdata(:,3);
	selected = userdata(:,4);
	peakanalysis = userdata(:,5);
else
	disp('No SSym XY data')
	return
end

if isempty(cmd) | strcmp(cmd,'get')
	if datfit == 2
		yapp = feval(appfcn,xapp,pin);
	end
% now verify that request is of good step, else interpolate
  if ~isempty(xin)
		dxin = abs(xin(1)-xin(2));
		dx = abs(x(1) - x(2));
		if abs(log(dxin/dx)) > 1e-3	% not full data
			nx = linspace(xapp(1),xapp(end),round(length(yapp)*dx/dxin));
			ny = interp1(xapp,yapp,nx,'nearest');
			nx = 1:length(ny);
			[dummy,centapp] = max(ny);	
			wapp = min(abs(centapp), abs(length(ny) - centapp));
			m1 = max(1,centapp - wapp);
			m2 = min(length(ny),centapp + wapp);
			ny = ny(m1:m2);
			yapp = yapp(:);
		end
	end
	yapp = yapp/sum(yapp); 
	return
end


% commands ----------------------------------------

if findstr(cmd,'hide')
	set(hsp_wapp,'Visible','off');
end

if findstr(cmd,'show')
	set(hsp_wapp,'Visible','on');
	figure(hsp_wapp);
end

if ~isempty(cmd)

	h=findobj(hsp_pars,'tag','hsp_ppeak');
	if ~isempty(h)
		peak = str2num(get(h,'string'));
	else
		disp('Can''t get Peak matrix')
		return
	end
	h=findobj(hsp_ctrl,'Tag','sp_xoap');
	if ~isempty(h)
		xoap = get(h,'userdata');
	else
		disp('Can''t get Xoapp')
		return
	end

	if datfit == 2
		[dummy, name2, pnames, pin] = feval(appfcn,xapp,[],1);	% identify
		pnames = strcat('App:',pnames);
%		pin = [];
%		name2 = appfcn;
	else
		pin = [];
		pnames = '';
		name2 = 'Data';
		delete(findobj('Tag','hsp_fitapp'));
	end
	name = [ name '+' name2 ];
	set(hsp_wapp,'Name',[ 'SSym : App Fcn with ' name2 ]);

end

if (~isempty(findstr(cmd,'fit')) | ~isempty(findstr(cmd,'guess'))) & (datfit == 2)

%	[dummy, name2, pnames, pin] = feval(appfcn,xapp,[],1);	% identify
%	pnames = strcat('App:',pnames);
	hmf_pars = findobj('Tag','mf_ParWindow');
	afmode = findobj(hsp_wapp,'tag','hsp_afauto');
	if isempty(afmode)
		afmode = 1;	% 1=auto, 2=manual
	else
		afmode = get(afmode,'Value');
	end

	fitrout = findobj('tag','mf_FitRoutineFile');
	if ~isempty(fitrout)
		fitrout=get(fitrout,'string');
	else
		fitrout = 'mf_flsqr';
	end
	i = find(x(peak(:,1)) == xoap);
	if isempty(i)
		disp('Can''t find App peak in matrix')
		return
	end
	peakapp = peak(i,:);
	peakapp(4) = peakapp(4) +xoap;
	if afmode == 1 | ~isempty(findstr(cmd,'guess'))
		disp('Guessing starting parameters')
		pin = mf_autogs(xapp,yapp,[],[],peakapp,pnames,appfcn);		% guess
		[ dummy, dummy, dummy, pin] = feval(appfcn,xapp,pin,3);
	elseif ~isempty(hmf_pars) & hmf_pars
		[pin, dp, fixed]=mf_rpars;
	else
		disp('No starting parameters !')
		return;
	end

	if afmode == 2
		newwin = get(hmf_pars,'Name');
		newwin = findstr(newwin,'SSym');
		if isempty(newwin) & ~isempty(hmf_pars) & hmf_pars
			[oldp, olddp] = mf_rpars;
			if ~isempty(oldp) & ~isempty(olddp)
				set(hsp_wapp,'userdata',[ oldp(:) olddp(:) ]);
			else
				set(hsp_wapp,'userdata',[]);
			end
		end
		delete(hmf_pars);
		hmf_pars = mf_pwin(name,pnames,[]);
		mf_upars(pin,[],[]);
		[pin, dp, fixed] = mf_rpars;
		set(hmf_pars,'Name',[ 'SSym : Parameters : ' name2 ]);
	else 
		fixed = [];
	end

	if ~isempty(findstr(cmd,'fit')) | (afmode == 1)

		[pin,dp] = feval(fitrout,xapp,yapp,errapp,pin,~fixed,appfcn,[]);		% fit
		fapp = feval(appfcn,xapp,pin);
		if ~isempty(hmf_pars) & hmf_pars & isempty(findobj('Tag','mf_ControlWindow'))
			mf_upars(pin,dp,[]);
		end
	end
	if ~isempty(pin) & ~isempty(find(pin))
		h = findobj(hsp_app,'Tag','hsp_fitapp');
		if ~isempty(h)
			delete(h);
		end
		wascurrent = get(0,'CurrentFigure');
		figure(hsp_ctrl);
		axes(hsp_app);
		hold on
		h = plot(xapp,feval(appfcn,xapp,pin),...
			'Color','red',...
			'erasemode','background',...
			'Tag','hsp_fitapp',...
			'LineStyle','--');
		hold off
% now fix parameters in MFit Params Window ------------------
		if ~isempty(dp)
			disp([ '* Apparatus function fit with ' appfcn ]);
			disp('---------------------------------------------')
			disp('  Parameter       Value      Uncertainty')

			for i=1:length(pin)
				disp(sprintf('%12s  %12.4e  %12.4e', pnames(i,:), pin(i), dp(i)))
			end
			disp('---------------------------------------------')
		else
			dp = 0*pin;
		end

		if ~isempty(hmf_pars) & hmf_pars & (~isempty(findstr(cmd,'fit')) | (afmode == 1))
			if ~isempty(findobj('Tag','mf_ControlWindow')) & ~isempty(findstr(get(hmf_pars,'Name'),'SSym'))
				mf_newfn;	% restore mfit parwin
				newwin = get(hsp_wapp,'userdata');
				if ~isempty(newwin)
					oldp = newwin(:,1);
					if size(newwin,2) > 1
						olddp = newwin(:,2);
					else 
						olddp = [];
					end
					newwin = [];
					set(hsp_app,'userdata',[]);
					mf_upars(oldp, olddp);
				end
			end
% now transfert values if needed -----------------------
			hpars = get(hmf_pars,'userdata');
			for i = 1:size(hpars,1)
				hpname = get(hpars(i,3),'String');
				for j = 1:size(pnames,1)
					if ~isempty(findstr(pnames(j,:),hpname))
						hpname = lower(hpname);
						n = [ findstr(hpname,'backg') findstr(hpname,'bkg') findstr(hpname,'const') findstr(hpname,'energ') findstr(hpname,'pos') findstr(hpname,'cent') findstr(hpname,'freq')];
						if ~isempty(n)	% tries to set background and centre to 0
							pin(j) = 0;
							fprintf(1,'Parameter %s is set to zero\n',pnames(j,:));
						end

						set(hpars(i,3),'value',1);
						set(hpars(i,1),'String',num2str(pin(j)));
						set(hpars(i,2),'String',num2str(dp(j)));
					end
				end
			end
		end
		figure(wascurrent)
	end
%	pnames = [];
%	pin = [];
end

if nargout == 0
	yapp = [ cmd ' done' ];
end


pin = pin(:);
