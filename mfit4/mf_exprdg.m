function [xexpr, yexpr,eexpr,xlab,ylab]=mf_exprdg(name)
%
% MFIT

if nargin == 0, name = ''; end
if isempty(name)
  name = 'MFIT : Direct X,Y,err';
end
%=========== Make column selection window... =======================

hmf_load=findobj('Tag','mf_exprdlg');
if ~isempty(hmf_load)
   figure(hmf_load);
   delete(gca);
   newflag=0;
else
   newflag=1;
%   mf_msg('Building Expr window');

%------- Create figure window -------------------------------------
   hmf_load=figure('Position',[400 400 300 151],...
      'Tag','mf_exprdlg',...
      'MenuBar','none',...
        'WindowStyle','modal',...
      'Name',name,...
      'Color',get(0,'DefaultUicontrolBackgroundColor'), ...
      'Resize','off',...
         'Visible','off',...
      'NumberTitle','off');

%--------  label ------------------------------------------
   uicontrol(hmf_load,...
      'Style','text',...
      'ForegroundColor',[1 1 1],...
      'Position',[10 132 300 16],...
      'String','Enter expressions for :');

%--------- edit labels -------------------------------------
   h(1)=uicontrol(hmf_load,...                       % x column label
      'Position',[5 106 50 22],...
      'String','x data');
   h(2)=uicontrol(hmf_load,...                       % y column label
      'Position',[5 80 50 22],...
      'String','y data');
   h(3)=uicontrol(hmf_load,...                       % error column label
      'Position',[5 54 50 22],...
      'String','y error');
   set(h,'Style','text','ForegroundColor',[1 1 1],'HorizontalAlignment','left');
   clear h;

%------ edit fields ------------------------------------------------

   h(1)=uicontrol(hmf_load,...                           % x column pop-up
      'Tag','mfexpr_xcol',...
      'Position',[60 106 200 25]);
   h(2)=uicontrol(hmf_load,...                           % y column pop-up
      'Tag','mfexpr_ycol',...
      'Position',[60 80 200 25]);
   h(3)=uicontrol(hmf_load,...                           % err column pop-up
      'Tag','mfexpr_ecol',...
      'Position',[60 54 200 25]);
   set(h,'Style','edit',...                             % Things common to all
      'BackgroundColor',[1 1 1],...
      'ForegroundColor',[0 0 0],...
      'HorizontalAlignment','right',...
    'Visible','on',...
    'Clipping','off');

%------- ok and cancel buttons -------------------------------------------
   uicontrol(hmf_load,...
      'Style','PushButton',...
      'String','Ok',...
      'Tag','mfexpr_ok',...
      'Position',[32 4 50 20]);
   uicontrol(hmf_load,...
      'Style','PushButton',...
      'String','Cancel',...
      'Tag','mfexpr_cancel',...
      'CallBack','delete(gcf)', ...
      'Position',[92 4 50 20]);

end
set(findobj('Tag','mf_exprdlg'),'Name',name);
set(findobj('Tag','mf_exprdlg'),'Visible','on');          % Switch window on


%============ Now wait until ok or cancel pressed =========================

drawnow;

hok = findobj('Tag','mfexpr_ok');
hcan= findobj('Tag','mfexpr_cancel');

if(0)
    waitforbuttonpress;
    while (hok & hcan & (gco ~= hok) & (gco ~= hcan))
        drawnow
        %   waitforbuttonpress;
    end
end


if(1)
    while (isempty(gco))
       drawnow
    end


%     while( ~isempty(gco) & (gco ~= hok) & (gco ~= hcan))
    while((hok ~= gco) & (hcan ~=gco ))
        drawnow
    end

    if(isempty(gco))
        axel1
    end

    my_gco = gco;
    delete(gco);

    uicontrol(hmf_load,...
       'Style','PushButton',...
       'String','Ok',...
       'Tag','mfexpr_ok',...
       'Position',[32 4 50 20]);

  end
%========= Extract relevant data columns ==================================

if (gco == hok)

%---------- Retrieve column numbers and labels from figure ---------------

   xexpr=get(findobj('Tag','mfexpr_xcol'),'String');          % menu positions
   yexpr=get(findobj('Tag','mfexpr_ycol'),'String');
   eexpr=get(findobj('Tag','mfexpr_ecol'),'String');

   xlab=xexpr(1:min(length(xexpr),10));
   ylab=yexpr(1:min(length(yexpr),10));

   set(findobj('Tag','mf_exprdlg'),'Visible','off');

else
   xexpr='';
   yexpr='';
   eexpr='';
   xlab='';
   ylab='';
   delete(findobj('Tag','mf_exprdlg'));
end

%mf_msg('Done');

