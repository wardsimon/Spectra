function mf_btprs
%
% MFIT function  mf_btprs
%	Process button press on data window
%	MZ 14.3.95
%

[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

%figure(hmf_data);
set(0,'CurrentFigure',hmf_data)

objtype=get(gco,'type');
presstype=get(gcf,'SelectionType');
% disp(['button press on: ' objtype '    type: ' presstype])


%====== Left mouse button press ========================================

if strcmp(presstype,'normal')

%------ if on text, edit it ---------------------
   if strcmp(objtype,'text')
   	mf_text(gco);

%------- if on axis, make zoom box ---------------
   elseif strcmp(objtype,'axes')
   	mf_rbbox('go');
   end

elseif strcmp(presstype,'alt') & strcmp(objtype,'axes')
	mf_rbbox('zoomout');


%====== Right mouse button press =======================================

elseif strcmp(presstype,'extend')

%------ if on axis, print coords of data point -------------------------
   if (strcmp(objtype,'axes') | strcmp(objtype,'line'))

      p=get(gca,'CurrentPoint');        % Work out which is closest point
    	data=get(hmf_data,'userdata');
   	x=data(:,1)';
   	y=data(:,2)';
      [d i]=min((x-p(1)).^2);                    % i is index of closest point

      txt=[' (' num2str(x(i)) ',' num2str(y(i)) ')'];
      axcol = get(findobj('tag','mf_AxesColor'),'string');
      if isempty(axcol), axcol = 'white'; end
      h=text(x(i),y(i),txt,...
           'Rotation',90,...
           'Fontsize',10,...
           'Tag','mf_point_coords',...
           'color',axcol,...
           'Fontweight','light');
      set(h,'userdata',h);
   end

%------ if on text, switch dragging on ------------------------
   if (strcmp(objtype,'text'))
      set(get(gco,'userdata'),'units','data','erasemode','xor');   % Set all erasemodes to xor
      set(gcf,'Pointer','Fleur');                                  % Cursor to fleur
      set(gcf,'WindowButtonMotionFcn','mf_drag');                  % Drag when moved
      set(gcf,'WindowButtonUpFcn',...                              % Stop when button released
         ['set(gcf,''WindowButtonMotionFcn'',''''),'...
          'set(gcf,''Pointer'',''arrow''),'...
          'set(get(gco,''userdata''),''erasemode'',''normal''),'...
          'set(gcf,''WindowButtonUpFcn'','''')']);
   end
end
