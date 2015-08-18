function [h]=mf_text(h, cmd)
% MFIT function mf_text(h)
%
% Complicated function for all text editing and adding in MFIT


%   Produce text editing dialog.
%     h is the handle of a text object to be edited.
%     (eg. axis and graph titles, user text)
%     MZ 29.11.94

[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;
%===== Process button presses ==================================
Str = '';
if nargin>1

%------ ok button pressed - do update ---------------------------
   if strcmp(cmd,'mf_text_ok')
    hh=get(gcf,'userdata');              % get handles
    h=get(hh(1),'userdata');             % handles of lines to update
    htxt=hh(2);                          % handle of edit window

      Str=get(htxt,'String');              % String to put
      mf_text(h(1),Str);                   % Handle of first element
      set(findobj('tag','mf_textdlg'),'visible','off');
      delete(h);

%------ font button pressed - select font -----------------------
   elseif strcmp(cmd,'mf_text_font')
    h=get(gcf,'userdata');
    uisetfont(h(1),'Select font:');
    hh=get(h(1),'userdata');
    for i=1:length(hh)
       set(hh(i),'FontName',get(h(1),'FontName'));
       set(hh(i),'FontSize',get(h(1),'FontSize'));
    end
      Str=get(h(2),'String');              % String to put
      mf_text(h(1),Str);                   % Handle of first element
      delete(hh);
      set(findobj('tag','mf_textdlg'),'visible','off');

%------ delete button pressed - delete object -------------------
   elseif strcmp(cmd,'mf_text_delete')
    h=get(gcf,'userdata');
  if ~isempty(h)
      hh=get(h(1),'userdata');
    if ~isempty(hh)
        delete(hh);
    end
  end
        set(findobj('tag','mf_textdlg'),'visible','off');

   elseif strcmp(cmd,'mf_text_deleteall')
  figure(hmf_data);
  hh=get(gca,'Children');
  for i=1:length(hh)
    h=hh(i);
    if (strcmp(get(h,'tag'),'mf_text_user') | strcmp(get(h,'tag'),'mf_point_coords'))
      delete(h);
    end
  end

%===== Add text requested =======================================
   elseif strcmp(cmd,'mf_text_add')

      nopars = 0;
      [hmf_ctrl, hmf_data, hmf_pars]=mf_figs;
      if (hmf_data==0 | hmf_pars==0)
        nopars = 1;
      end

%------Generate string matrix of parameter names and values------
      if (nopars == 0)
        h=get(hmf_pars,'userdata');
        pnames=[];
        pvals=[];
        perr=[];
        for i=1:size(h,1)
           pnames=str2mat(pnames,[get(h(i,3),'String') '  ']);
           pvals =str2mat(pvals ,[get(h(i,1),'String') ' +/- ']);
           perr  =str2mat(perr  ,get(h(i,2),'String'));
        end
        [r2,rv] = mf_stats;
        time = clock;
        str=[pnames pvals perr];
        str(1,:)=[];
        str=str2mat([ 'FitFcn: ' get(findobj('Tag','mf_FitFuncName'),'String')...
                    ' (' get(findobj('Tag','mf_FitFuncFile'),'String') ')' ],...
                    str,['Chi^2       ' get(findobj('tag','ChiSq'),'string')]);
        str=str2mat(str,sprintf('Res. Var.  %.3f',rv), sprintf('Cor. Coef. %.4f',r2), sprintf('Date :     %s %i:%i\n', date, time(4), time(5)));
        hmf_convlv = findobj('Tag', 'mf_convlv');
       if ~isempty(hmf_convlv)
        UserData = get(hmf_convlv, 'UserData');
        str = str2mat(str, char(UserData.input.title));
       end
      else
        str='';
      end
      toadd = get(findobj('Tag','mf_DataDir'),'String');
      i = findstr(toadd,'\');
      if ~isempty(i)
        toadd(i) = '/'; % avoid TeX string in path
      end
      str=str2mat(' ',[ 'Data  : ' get(findobj('Tag','mf_DataFile'),'String') ],...
                      [ '   in : ' toadd ],...
                  str);
      figure(hmf_data);
      hax=get(hmf_data,'currentaxes');
      h=text(0,0,'','units','normalized',...
            'position',[0 1 0],...
            'Fontsize',10,...
            'FontName','courier',...
            'Tag','mf_text_user');

      h=mf_text(h,str);
      set(get(h,'userdata'),'visible','off');
      mf_text(h);

   elseif h

%------ place text requested ------------------------------------
      Str=cmd;                                  % Matlab string to place
      hax=get(h,'parent');                      % Axes to put it on
      set(h,'units','points');
      pos=get(h,'position');                    % Get origin of text in points
      Size=get(h,'fontsize');                   % Get fontsize
      axes(hax);
      if (size(Str,1)  == 0)
  Str = ' ';
      end
      axcolor = get(findobj('tag','mf_AxesColor'),'string');
      if isempty(axcolor) axcolor = 'white'; end
      for i=1:size(Str,1)
         hh(i)=text(0,0,Str(i,:),...            % Place line by line, (vec. doesn't work)
               'FontSize',Size,...              % retaining template properties'
               'color',axcolor,...
               'Tag',get(h,'tag'),...
               'Units','points',...
               'FontName',get(h,'FontName'),...
               'Rotation',get(h,'Rotation'),...
               'HorizontalAlignment',get(h,'HorizontalAlignment'),...
               'VerticalAlignment',get(h,'VerticalAlignment'));

         th=pi/180*get(h,'rotation');           % Work out position for new line
         R=[cos(th) -sin(th); sin(th) cos(th)]; % Rotation matrix
         p=[0; Size*(i-1)];                     % Offset in y
         pp=R*p;                                % New offset vector
         set(hh(i),'position',...
             [pos(1)-pp(1) pos(2)-pp(2) 0]);    % Set position = orig + offset
      end
      set(hh,'userdata',hh);      % Set all userdata's to hold
      h=hh(1);          % handles of other objects
   end


%===== else make edit text dialog =================================

else

   if isempty(findobj('Tag','mf_textdlg'))

    Title='MFIT: Edit text:';
    boxh=100;
    Min=0; Max=2;
    TextBoxHeight = get(findobj('tag','mf_TextBoxHeight'),'string');
  if isempty(TextBoxHeight)
    TextBoxHeight=18;
  else
    TextBoxHeight = str2num(TextBoxHeight);
  end

    htext=figure('Position',[5 320 400 TextBoxHeight+6+boxh],...
               'Name',Title,...
               'NumberTitle','off',...
               'MenuBar','none',...
               'Color',[0.4 0.4 0.4],...
               'BackingStore','off',...
               'Resize','off',...
               'Tag','mf_textdlg',...
               'Parent',0,...
               'visible','off');

    htxt=uicontrol(htext,...
               'Style','edit',...
               'String',Str,...
               'Position',[2 TextBoxHeight+4 396 boxh],...
               'BackgroundColor',[1 1 1],...
               'ForegroundColor',[0 0 0],...
               'Min',Min,...
               'Max',Max,...
               'HorizontalAlignment','Left');

%--- Font button ------------------------------------------
      uicontrol(htext,...
               'Style','push',...
               'String','Font...',...
               'Position',[20 2 70 TextBoxHeight],...
               'Callback','mf_text(0,''mf_text_font'');');
%--- Delete button ------------------------------------------
      uicontrol(htext,...
               'Style','push',...
               'String','Delete',...
               'Position',[110 2 70 TextBoxHeight],...
               'Callback','mf_text(0,''mf_text_delete'');');

%--- Ok button ------------------------------------------
      uicontrol(htext,...
               'Style','push',...
               'String','Ok',...
               'Position',[220 2 70 TextBoxHeight],...
               'Callback','mf_text(0,''mf_text_ok'');');

%--- Cancel button ------------------------------------------
      uicontrol(htext,...
               'Style','push',...
               'String','Cancel',...
               'Position',[310 2  70 TextBoxHeight],...
               'Callback','set(gcf,''visible'',''off'');');

      set(htext,'userdata',[h htxt]);
   end

%===== Put text into dialog =======================================

   htext=findobj('Tag','mf_textdlg');         % Get edit figure handle
   htext = htext(1);
   set(htext,'Visible','on');                 % Switch on

%----- Make text from handle passed ------------------------------
   hh=get(h,'userdata');                      % Handles of all linked lines
   if isempty(hh)
      set(h,'userdata',h);
      hh=h;
   end
   Str=get(hh(1),'String');
   for i=2:length(hh)
      Str=str2mat(Str,get(hh(i),'String'));   % Make matlab string of lines
   end

%----- Set edit window string ------------------------------------
   hh=get(htext,'userdata');
   if length(hh) > 1
    htxt=hh(2);                                % Handle of edit window
    set(htxt,'String',Str);
    set(htext,'userdata',[h htxt]);
   end
end
