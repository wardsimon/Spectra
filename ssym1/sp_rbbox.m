function sp_rbbox(m)
%
% SSym function sp_rbbox(m)
% 		Handle zooming and data selection and object selection
% 		EF 04.97 from MZ 29.11.94

% m='zoomin'		Zoom in
% m='zoomout'		Zoom out
% m='zoomreset'		Zoom reset
% m='select'		Select data/peak
% m='selectpk'		Select peak
% m='selectbg'		Select background
% m='selectapp'		Select app zone
% m='deselect'		Deselect data
% m='selectall'  	Select all data
% m='deselectall'	Deselect all data

index=[]; cmd = '';
markersize = 6;
if ~exist('m')
	m = '';
end

hsp_ctrl = findobj('tag','hsp_ctrl');
userdata = get(hsp_ctrl,'userdata');
hsp_plot = userdata(1);
hsp_app  = userdata(2);

figure (hsp_ctrl);
hax = hsp_plot;
x = []; y=[]; err=[]; selected = []; pka = [];
presstype=get(gcf,'SelectionType');

cmd=get(findobj('tag','sp_rbbox_action'),'userdata');
if findstr(cmd,'zoom')
	set(findobj('Tag','hsp_mvzoom'),'Checked','on');
else
	set(findobj('Tag','hsp_mvzoom'),'Checked','off');
end
if strcmp(cmd,'select')
	set(findobj('Tag','hsp_mvsel'),'Checked','on');
else
	set(findobj('Tag','hsp_mvsel'),'Checked','off');
end
if strcmp(cmd,'selectpk')
	set(findobj('Tag','hsp_mvsep'),'Checked','on');
else
	set(findobj('Tag','hsp_mvsep'),'Checked','off');
end
if strcmp(cmd,'selectbg')
	set(findobj('Tag','hsp_mvseb'),'Checked','on');
else
	set(findobj('Tag','hsp_mvseb'),'Checked','off');
end
if strcmp(cmd,'selectapp')
	set(findobj('Tag','hsp_mvsea'),'Checked','on');
else
	set(findobj('Tag','hsp_mvsea'),'Checked','off');
end
if strcmp(cmd,'deselect')
	set(findobj('Tag','hsp_mvusel'),'Checked','on');
else
	set(findobj('Tag','hsp_mvusel'),'Checked','off');
end

if strcmp(m,'go') & strcmp(presstype,'alt')
	sp_rbbox('zoomout');
	return;
end

if strcmp(m,'go') & strcmp(presstype,'extend')
	fprintf(1,'The current mouse position is : ');
	[x,y] = ginput(1);
	fprintf(1,'[ %g , %g ]\n',x,y);
	return
end

%----------- Extract data from figure ---------------
data=get(hsp_plot,'userdata');
if ~isempty(data)
	x=data(:,1);
	y=data(:,2);
	err=data(:,3);
	index=data(:,4);
	pka = data(:,5);
end
		
%-------- Store current action for rbbox in rbbox action -----------
if strcmp(m,'zoomin') | strcmp('select',m) | strcmp('deselect',m) ...
	| strcmp('selectbg',m) | strcmp('selectpk',m) | strcmp('selectapp',m)
	h=findobj('tag','sp_rbbox_action');
	set(h,'userdata',m);
	return;
end

%----------- Zoom out -----------------------------------
if strcmp(m,'zoomout') 
	h=findobj('tag','sp_zoom_list');	% get old zooms
	oldzooms=get(h,'userdata');
	n=size(oldzooms,1);
	if n>=1
		pos=oldzooms(n,:);
		zooms=oldzooms(1:n-1,:);
		set(h,'userdata',zooms);
		set(hax,'XLim',[pos(1) pos(2)]);		% Set new position
		set(hax,'YLim',[pos(3) pos(4)]);
	end

%----------- Zoom reset -----------------------------------
elseif strcmp(m,'zoomreset') 
	h=findobj('tag','sp_zoom_list');					
	set(h,'userdata',[]);
	axes(hsp_plot)
	axis auto
	drawnow
	tmp = axis;
	tmp(1) = min(x); tmp(2) = max(x);
	axis(tmp);
	sp_rbbox('zoomin');

%--------- Select/Deselect all ------------------
elseif strcmp('selectall',m)
	index=ones(size(index));
	pka = abs(pka);
	set(hsp_plot,'userdata',[x y err pka index]);
elseif strcmp('deselectall',m)
	index=zeros(size(index));
	pka = -abs(pka);
	set(hsp_plot,'userdata',[x y err pka index]);

%--------- Called by button down event... -----------
elseif strcmp('go',m)

%----------- Get box coordinates ----------------------
	pt1=get(hax,'CurrentPoint');			% save current point
	p=get(gcf,'CurrentPoint');
	rbbox([p 0 0],p);							% call rbbox
	pt2=get(hax,'CurrentPoint');			% get new point
	xmin=min([pt1(1,1) pt2(1,1)]);
	xmax=max([pt1(1,1) pt2(1,1)]);
	ymin=min([pt1(1,2) pt2(1,2)]);
	ymax=max([pt1(1,2) pt2(1,2)]);

%----------- Retrieve stored command -------------------
	cmd=get(findobj('tag','sp_rbbox_action'),'userdata');
	if (xmin==xmax | ymin==ymax) return; end;

%----------- Zoom in -----------------------------------
	if strcmp(cmd,'zoomin') 
		h=findobj('tag','sp_zoom_list');	% Store old position
		oldzooms=get(h,'userdata');
		curax = [get(hax,'XLim') get(hax,'YLim')];
		xmin = min(curax(2),max(xmin,curax(1)));
		xmax = max(curax(1),min(xmax,curax(2)));
		ymin = min(curax(4),max(ymin,curax(3)));
		ymax = max(curax(3),min(ymax,curax(4)));
		if (xmin ~= xmax) & (ymin ~= ymax)
			set(hax,'XLim',[xmin xmax]);		% Set new position
			set(hax,'YLim',[ymin ymax]);
			zooms=[oldzooms; curax ];
			n=size(zooms,1);
			if (n>10) zooms=zooms(n-10:n,:); end;
			set(h,'userdata',zooms);
		end
	end

	if (strcmp(cmd,'select') | strcmp(cmd,'deselect'))
		sel=(x>xmin & x<xmax & y>ymin & y<ymax);
		if strcmp(cmd,'select')
			index=index | sel;
			sel = find(sel);
			pka(sel) = abs(pka(sel));
		else
			index=index & ~sel;
			sel = find(sel);
			pka(sel) = -abs(pka(sel));
		end
	end
	if strcmp(cmd,'selectpk')
		sel=find(x>xmin & x<xmax & y>ymin & y<ymax);
		pka(sel) = 2;
	end
	if strcmp(cmd,'selectbg')
		sel=find(x>xmin & x<xmax & y>ymin & y<ymax);
		pka(sel) = 3;
	end
	if strcmp(cmd,'selectapp')
		sel=find(x>xmin & x<xmax & y>ymin & y<ymax);
		pka(sel) = 1;
		hsp_pars = findobj('tag','hsp_pars');
		h=findobj(hsp_pars,'tag','hsp_pshut');
		if ~isempty(h)
			shutterzones = str2num(str2mat('[ ',get(h,'string'),' ]'));
			shutterzones = [ shutterzones ; xmin xmax ];
			set(h,'String',num2str(shutterzones,'%.8g '));
		end
	end
	set(hsp_plot,'userdata',[x y err index pka]);
end

curax = [get(hsp_plot,'XLim') get(hsp_plot,'YLim')];
i = find(x>curax(1) & x<curax(2) & y>curax(3) & y<curax(4));
if ~isempty(i)
	markersize2 = min(20,markersize*length(x)/length(i));
else
	markersize2 = markersize;
end
h = findobj(hsp_plot,'Tag','hsp_plsel');
set(h,'MarkerSize',markersize2);

if (findstr(cmd,'select') & strcmp(m,'go')) | findstr(m,'select')
	set(hsp_plot,'userdata',[x y err index pka ]);
	ssym;
end
	
