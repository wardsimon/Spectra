function mf_rbbox(m)
%
% MFIT function mf_rbbox(m)
%     Handle zooming and data selection and object selection
%     MZ 29.11.94
%
% m='zoomin'    Zoom in
% m='zoomout'   Zoom out
% m='zoomreset'   Zoom reset
% m='select'    Select data
% m='deselect'    Deselect data
% m='selectall'   Select all data
% m='deselectall' Deselect all data

index=[]; cmd = '';
if nargin == 0
  m = '';
end
[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

figure (hmf_data);

hax=get(hmf_data,'CurrentAxes');

%-------- Store current action for rbbox in rbbox action -----------
if (strcmp(m,'zoomin') | strcmp('select',m) | strcmp('deselect',m))
  h=findobj('tag','mf_rbbox_action');
  if (length(h) > 1)
    for i=2:length(h)
      delete(h(i));
    end
  end
  h = h(1);
  set(h,'userdata',m);
  return;
end

%----------- Extract data from figure ---------------

data=get(hmf_data,'userdata');
if ~isempty(data)
  x=data(:,1);
  y=data(:,2);
  err=data(:,3);
  index=data(:,4);
end

%----------- Zoom out -----------------------------------
if strcmp(m,'zoomout')
  h=findobj('tag','mf_zoom_list');  % get old zooms
  oldzooms=get(h,'userdata');
  if iscell(oldzooms), oldzooms=[ min(x) max(x) min(y) max(y) ]; end
  n=size(oldzooms,1);
  if n>=1
    pos=oldzooms(n,:);
    zooms=oldzooms(1:n-1,:);
    set(h,'userdata',zooms);
    set(hax,'XLim',[pos(1) pos(2)]);    % Set new position
    set(hax,'YLim',[pos(3) pos(4)]);
  else
    mf_rbbox('zoomreset')
  end

%----------- Zoom reset -----------------------------------
elseif strcmp(m,'zoomreset')
  h=findobj('tag','mf_zoom_list');  % Store old position
  if (length(h) > 1)
    for i=2:length(h)
      delete(h(i));
    end
  end
  h = h(1);
  select = find(index);
  oldzooms=[ min(x(select)) max(x(select)) min(y(select)) max(y(select)) ];
  pos = oldzooms;
  set(hax,'XLim',[pos(1) pos(2)]);    % Set new position
  set(hax,'YLim',[pos(3) pos(4)]);
  mf_rbbox('zoomin');

%--------- Select/Deselect all ------------------
elseif strcmp('selectall',m)
  index=ones(size(index));
elseif strcmp('deselectall',m)
  index=zeros(size(index));

%--------- Called by button down event... -----------
elseif strcmp('go',m)

%----------- Get box coordinates ----------------------
  pt1=get(hax,'CurrentPoint');      % save current point
  p=get(gcf,'CurrentPoint');
  rbbox([p 0 0],p);             % call rbbox
  pt2=get(hax,'CurrentPoint');      % get new point
  xmin=min([pt1(1,1) pt2(1,1)]);
  xmax=max([pt1(1,1) pt2(1,1)]);
  ymin=min([pt1(1,2) pt2(1,2)]);
  ymax=max([pt1(1,2) pt2(1,2)]);

%----------- Retrieve stored command -------------------
  cmd=get(findobj('tag','mf_rbbox_action'),'userdata');
  if iscell(cmd), cmd = 'zoomin'; end
  if (xmin==xmax | ymin==ymax) return; end;

%----------- Zoom in -----------------------------------
  if strcmp(cmd,'zoomin')
      h=findobj('tag','mf_zoom_list');          % Store old position
    oldzooms =get(h,'userdata');
    if isempty(oldzooms)
            oldzooms=[ min(x) max(x) min(y) max(y) ] ;
        end
%     oldzooms
%          get(hax,'XLim')
%          get(hax,'YLim')
        zooms=[oldzooms; [get(hax,'XLim') get(hax,'YLim')]];
    n=size(zooms,1);
    if (n>10) zooms=zooms(n-10:n,:); end;
    set(h,'userdata',zooms);
    set(hax,'XLim',[xmin xmax]);    % Set new position
    set(hax,'YLim',[ymin ymax]);
  end

  if (strcmp(cmd,'select') | strcmp(cmd,'deselect'))
    sel=(x>xmin & x<xmax & y>ymin & y<ymax);
    if strcmp(cmd,'select')
      index=index | sel;
    else
      index=index & ~sel;
    end
  end
end

cmd=get(findobj('tag','mf_rbbox_action'),'userdata');
if iscell(cmd), cmd = 'zoomin'; end
if findstr(cmd,'zoom')
  set(findobj('Tag','hmf_mvzoom'),'Checked','on');
else
  set(findobj('Tag','hmf_mvzoom'),'Checked','off');
end
if strcmp(cmd,'select')
  set(findobj('Tag','hmf_mvsel'),'Checked','on');
else
  set(findobj('Tag','hmf_mvsel'),'Checked','off');
end
if strcmp(cmd,'deselect')
  set(findobj('Tag','hmf_mvusel'),'Checked','on');
else
  set(findobj('Tag','hmf_mvusel'),'Checked','off');
end

markersize = str2num(get(findobj('Tag','mf_MarkerSize'),'String'));
if isempty(markersize)
  markersize = 10;
end

limits=[get(hax,'Xlim') get(hax,'Ylim')];
i = find(x>limits(1) & x<limits(2) & y>limits(3) & y<limits(4));
if ~isempty(i)
  t = log(length(i)/100); t = max(t,.1);
  markersize = round(markersize / t);
  markersize = max(8,min(20, markersize));
end

set(findobj('Tag','mf_selected'),'MarkerSize',markersize);
set(findobj('Tag','mf_alldat'),'MarkerSize',ceil(markersize/2));

if findstr('select',[ cmd m ])
  set(hmf_data,'userdata',[x y err index]);
  mf_uplot('sel');
end
%drawnow



