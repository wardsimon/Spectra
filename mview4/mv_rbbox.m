function mv_rbbox(m)
%
% MFIT  function mv_rbbox(m)
%   Handle zooming and data selection and object selection
%   MZ 29.11.94

% m='zoomin'    Zoom in
% m='zoomout'   Zoom out
% m='zoomreset' Zoom reset
% m='select'    Select data
% m='deselect'    Deselect data
% m='selectall'   Select all data
% m='deselectall' Deselect all data

index=[]; cmd = '';
[hmv_ctrl, hmv_data]=mv_figs;

hax=get(hmv_data,'CurrentAxes');

%-------- Store current action for rbbox ------------
if (strcmp(m,'zoomin') | strcmp('select',m) | strcmp('deselect',m))
  h=get(hax,'XLabel');
  set(h,'userdata',m);
  return;
end

%----------- Extract data from figure ---------------
if (strcmp('go',m) | strcmp('selectall',m) | strcmp('deselectall',m))
  data=get(hmv_data,'userdata');
  x=data(:,1);
  y=data(:,2);
  err=data(:,3);
  index=data(:,4);
end

%----------- Zoom out -----------------------------------
if strcmp(m,'zoomout')
  h=get(hax,'ZLab');          % Store old position
  oldzooms=get(h,'userdata');
  n=size(oldzooms,1);
  if n>=1
    pos=oldzooms(n,:);
    zooms=oldzooms(1:n-1,:);
    set(h,'userdata',zooms);
    set(hax,'XLim',[pos(1) pos(2)]);    % Set new position
    set(hax,'YLim',[pos(3) pos(4)]);
  end

%----------- Zoom reset -----------------------------------
elseif strcmp(m,'zoomreset')
  h=get(hax,'ZLab');          % Store old position
  oldzooms=get(h,'userdata');
  if ~isempty(oldzooms)
    pos=oldzooms(1,:);
    set(h,'userdata',[]);
    set(hax,'XLim',[pos(1) pos(2)]);    % Set new position
    set(hax,'YLim',[pos(3) pos(4)]);
  end

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
  rbbox([p 0 0],p);       % call rbbox
  pt2=get(hax,'CurrentPoint');      % get new point
  xmin=min([pt1(1,1) pt2(1,1)]);
  xmax=max([pt1(1,1) pt2(1,1)]);
  ymin=min([pt1(1,2) pt2(1,2)]);
  ymax=max([pt1(1,2) pt2(1,2)]);

%----------- Retrieve stored command -------------------
  cmd=get(get(hax,'Xlab'),'userdata');

  if (xmin==xmax | ymin==ymax) return; end;

%----------- Zoom in -----------------------------------
  if strcmp(cmd,'zoomin')
    h=get(hax,'ZLab');      % Store old position
    oldzooms=get(h,'userdata');
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

if ~isempty(findstr('select',[cmd m]))
  set(hmv_data,'userdata',[x y err index]);
        mv_uplot('sel');
else
   mv_uplot('all');
end




