function [xcolo, ycolo, ecolo, mcolo, xlab, ylab, scano, normo,keep,flipo ]=mf_coldg(collist, ncol,xcol,ycol,ecol,mcol,scan,normi,Fstate, coldg_title, file_index)
%
% MFIT column selector
%
% MZ 5.6.95 EF 15.09.97 (update) ARW 11.98 (polarization)
% ecol = n+1 for sqrt(y), n+2 for none
% mcol = n+1 for none
% returns also x and y labels, optional fileds, and keep = 1 if window
% musn't be hidden.

if nargin < 11, file_index=[]; end
if nargin < 10, coldg_title=[]; end
if nargin < 9 , Fstate=[]; end
if nargin < 8 , normi=[]; end
if nargin < 7 , scan=[]; end
if nargin < 6 , ecol=[]; end
if nargin < 9 , ycol=[]; end
if nargin < 5 , xcol=[]; end
if nargin < 4 , ncol=[]; end

%======= Generate menu list for x,y, err column pop ups =================
%
% This section generates the list of column names that will be offered to
% the user in the column pop-up menus. First see if the last line of the
% header in the data file might be column names (right number of words)...

   xcolo =[];
   ycolo=[];
   mcolo=[];
   ecolo=[];
   scano='';
   normo=[];
   flipo=[];
   keep = 0;

%-----Try to be clever and generate list from last line of header ----------
i=find(collist==9);                                   % replace tabs by spaces
collist(i)=32*ones(size(i));

[nncols,j]=size(collist);

%collist=abs(deblank(fliplr(deblank(fliplr(setstr(collist))))));    % strip leading and trailing spaces
%i=find(filter([1 1],2,collist==' ')==1);               % remove multiple spaces
%collist(i)=[];
%j=(collist==' ');
%nncols=sum(j)+1;                                      % work out number of columns

%----- Right number of columns? No? then just numbers...------------------
if nncols==ncol                                       % Generate number list
  ndatastr = collist(1,:);
  for i=2:ncol          % set a header line
    ndatastr = [ ndatastr '|' collist(i,:) ];
  end
  collist = ndatastr;
else
   collist=['column 1' sprintf('|column %2d',2:ncol)];
end

%if isempty(xcol) xcol = 1; end
%if isempty(ycol) ycol = 2; end
%if isempty(ecol) ecol = ncol+1; end
%if isempty(mcol) mcol = ncol+1; end

%=========== Increase size of window if polarisation is necessary
plist='???';

if ~isempty(Fstate)
  Fnum = 50;
  if min(size(Fstate)) > 1
    plist = sprintf('F1=%d,F2=%d',Fstate(1,1),Fstate(1,2));
    for i=2:size(Fstate(:,1))
      plist=[plist, sprintf('|F1=%d,F2=%d',Fstate(i,1),Fstate(i,2))];
    end
  else
    plist = sprintf('PAL=%d',Fstate(1));
    if length(Fstate) > 1
      plist = [ plist sprintf('|PAL=%d',Fstate(2:end)) ];
    end
  end
else
  Fnum = 0;
  Fstate=[0 0];
end

if isempty(coldg_title), coldg_title='Column Selector'; end

%=========== Make column selection window... =======================

hmf_load=findobj('Tag','mf_coldlg');
if ~isempty(hmf_load)
   figure(hmf_load);
   delete(gca);
   newflag=0;

   set(hmf_load, 'Name', coldg_title);
   set(hmf_load,'Position',[800 200 174 207 + Fnum]);   %--- Reset window size and button positions
   ypos=[106 80 54 28 134+20 162+20]+Fnum;
   hx=findobj('Tag','mfload_xcol');       set(hx,'Position',[60 ypos(1) 110 20]);
   hy=findobj('Tag','mfload_ycol');       set(hy,'Position',[60 ypos(2) 110 20]);
   he=findobj('Tag','mfload_ecol');       set(he,'Position',[60 ypos(3) 110 20]);
   hm=findobj('Tag','mfload_mcol');       set(hm,'Position',[60 ypos(4) 110 20]);
   hp=findobj('Tag','mfload_pol');        set(hp,'Position',[60 ypos(4)-Fnum/2 110 20],'Visible','on');
   hs=findobj('Tag','mfload_scan');       set(hs,'Position',[85 ypos(5) 80 20]);
   hn=findobj('Tag','mfload_norm');       set(hn,'Position',[85 ypos(6) 80 20]);
   hf1=findobj('Tag','mfload_F1');        set(hf1,'Position',[50 ypos(4)-Fnum 22 22],'Visible','on');
   hf2=findobj('Tag','mfload_F2');        set(hf2,'Position',[125 ypos(4)-Fnum 22 22],'Visible','on');

   hcl=findobj('Tag','mfload_clabel');    set(hcl,'Position',[90 132 50 16 + Fnum]);
   hxl=findobj('Tag','mfload_xlabel');    set(hxl,'Position',[5 ypos(1) 50 22]);
   hyl=findobj('Tag','mfload_ylabel');    set(hyl,'Position',[5 ypos(2) 50 22]);
   hel=findobj('Tag','mfload_elabel');    set(hel,'Position',[5 ypos(3) 50 22]);
   hml=findobj('Tag','mfload_mlabel');    set(hml,'Position',[5 ypos(4) 50 22]);
   hsl=findobj('Tag','mfload_spolabel');  set(hsl,'Position',[5 ypos(5) 80 16]);
   hynl=findobj('Tag','mfload_ynlabel');  set(hynl,'Position',[15 ypos(6) 50 16]);
   hpl=findobj('Tag','mfload_plabel');    set(hpl,'Position',[5 ypos(4)-Fnum/2 50 22],'Visible','on');
   hf1l=findobj('Tag','mfload_f1label');  set(hf1l,'Position',[25 ypos(4)-Fnum 22 22],'Visible','on');
   hf2l=findobj('Tag','mfload_f2label');  set(hf2l,'Position',[100 ypos(4)-Fnum 22 22],'Visible','on');

   if Fnum == 0
     set(hp,'Visible','off');
     set(hf1,'Visible','off');
     set(hf2,'Visible','off');
     set(hpl,'Visible','off');
     set(hf1l,'Visible','off');
     set(hf2l,'Visible','off');
   end

   % when dialog already exists, retain previous user choice if column names didn't change
   tmp1 = cellstr(get(hy,'String'));
   tmp1 = strcat(tmp1{:});
   tmp2 = strrep(strcat(collist), '|','');
   tmp2 = strrep(tmp2, ' ','');
   if strcmp(tmp1, tmp2)
     xcol = [];
     ycol = [];
     ecol = [];
     mcol = [];
     scan = [];
     normi= [];
  else
    file_index = [];
  end
else
   newflag=1;

%------- Create figure window -------------------------------------
   hmf_load=figure('Position',[800 200 174 207+Fnum],...
       'WindowStyle','modal',...
       'Tag','mf_coldlg',...
       'MenuBar','none',...
       'Name',coldg_title,...
       'Color',get(0,'DefaultUicontrolBackgroundColor'), ...
       'Resize','off',...
           'Visible','on',...
       'NumberTitle','off');

      %truc = get(hmf_load);

%-------- 'Column' label ------------------------------------------
   uicontrol(hmf_load,...
       'Style','text',...
       'ForegroundColor',[1 1 1],...
       'Position',[90 132 50 16+Fnum],...
       'String','Column', ...
       'ToolTipString','Please, assign the plot axes to data columns.');

ypos=[106 80 54 28 134+20 162+20]+Fnum;

%--------- column name labels -------------------------------------
   NL = sprintf('\n');
   h(1)=uicontrol(hmf_load,...                       % x column label
       'Position',[5 ypos(1) 50 22],...
    'Tag','mfload_xlabel',...
       'String','X data', ...
       'ToolTipString','Select the values to use for X axis');
   h(2)=uicontrol(hmf_load,...                       % y column label
       'Position',[5 ypos(2) 50 22],...
    'Tag','mfload_ylabel',...
       'String','Y data', ...
       'ToolTipString','Select the values to use for Y axis (signal)');
   h(3)=uicontrol(hmf_load,...                       % error column label
       'Position',[5 ypos(3) 50 22],...
    'Tag','mfload_elabel',...
       'String','Y error', ...
       'ToolTipString',['Select the values to use for Y errors' NL '(e.g. sqrt(y) or none)']);
   h(4)=uicontrol(hmf_load,...                       % monitor column label
       'Position',[5 ypos(4) 50 22],...
    'Tag','mfload_mlabel',...
       'String','Monitor', ...
       'ToolTipString',['Select the values to use for the Monitor' NL '(divides the Y values)']);
   h(5)=uicontrol(hmf_load,...                       % Polarisation state label
       'Position',[5 ypos(4)-Fnum/2 50 22],...
    'Tag','mfload_plabel',...
       'String','Pol state', ...
       'ToolTipString',[ 'Select the Polarisation State to use ' NL '(some part of Y values)' ]);

   set(h,'Style','text','ForegroundColor',[1 1 1],'HorizontalAlignment','left');
   if Fnum == 0
      set( h(5),'Visible','off');       % Get rid of flipper labels if they're not used
   end;

   clear h;

%------ Make pop-ups menus ------------------------------------------------

   h(1)=uicontrol(hmf_load,...                           % x column pop-up
       'Tag','mfload_xcol',...
       'Position',[60 ypos(1) 110 20],...
    'Value',ncol+1);
   h(2)=uicontrol(hmf_load,...                           % y column pop-up
       'Tag','mfload_ycol',...
    'Value',2,...
       'Position',[60 ypos(2) 110 20]);
   h(3)=uicontrol(hmf_load,...                           % err column pop-up
       'Tag','mfload_ecol',...
    'Value',ncol+1,...
       'Position',[60 ypos(3) 110 20]);
   h(4)=uicontrol(hmf_load,...                           % monitor column pop-up
       'Tag','mfload_mcol',...
    'Value',ncol+1,...
       'Position',[60 ypos(4) 110 20]);
   h(5)=uicontrol(hmf_load,...                           % Polarisation pop-up
    'Tag','mfload_pol',...
    'Value',1,...
    'Position',[60 ypos(4)-Fnum/2 110 20]);

   set(h,'Style','popup',...                             % Things common to all
       'BackgroundColor',[1 1 1],...
       'ForegroundColor',[0 0 0],...
       'HorizontalAlignment','right',...
    'Visible','on',...
    'Clipping','off');

%------- optional fields : nscan norm ------------------------------------

   uicontrol(hmf_load,...
          'Style','text',...
          'ForegroundColor',[1 1 1],...
          'Position',[5 ypos(5) 80 16],...
          'HorizontalAlignment','left',...
    'Tag','mfload_spolabel',...
          'String','scan/part/opt');

   h(6)=uicontrol(hmf_load,...
             'Style','edit',...
       'Tag','mfload_scan',...
       'String','1',...
       'Position',[85 ypos(5) 80 20],...
       'BackgroundColor',[1 1 1],...
       'ForegroundColor',[0 0 0],...
       'HorizontalAlignment','Left', ...
       'ToolTipString',[ 'Enter here a Scan sub-part or an option ' NL ...
       '(e.g. "1,setpar" to set Rescal parameters from ILL TAS files)' ]);

   uicontrol(hmf_load,...
          'Style','text',...
          'ForegroundColor',[1 1 1],...
          'Position',[15 ypos(6) 50 16],...
          'HorizontalAlignment','left',...
      'Tag','mfload_ynlabel',...
          'String','ynorm');

   h(7)=uicontrol(hmf_load,...
             'Style','edit',...
       'Tag','mfload_norm',...
       'String','1',...
       'Position',[85 ypos(6) 80 20],...
       'BackgroundColor',[1 1 1],...
       'ForegroundColor',[0 0 0],...
       'HorizontalAlignment','Left', ...
       'ToolTipString',['Enter here the normalisation to operate ' NL '(usually 1, or "1/sum(y)")']);

%------- Information on flipper states -----------------------------------


   if Fnum == 0
  set( h(5),'Visible','off');       % Get rid of polarisation pop-up if not used
   end;

%------- ok and cancel buttons -------------------------------------------
   uicontrol(hmf_load,...
       'Style','PushButton',...
       'String','Ok',...
       'Tag','mfload_ok',...
       'Position',[32 4 50 20], ...
       'ToolTipString',['Click here to import file ' NL ...
       'and save your settings' NL 'as long as the file type does not change']);
   uicontrol(hmf_load,...
       'Style','PushButton',...
       'String','Cancel',...
       'Tag','mfload_cancel',...
       'Position',[92 4 50 20], ...
       'ToolTipString','Click here to cancel the importation');

end

%------- Set menu lists and choices ----------------------------------------


hx=findobj('Tag','mfload_xcol');                         % Get current menu choices
hy=findobj('Tag','mfload_ycol');
he=findobj('Tag','mfload_ecol');
hm=findobj('Tag','mfload_mcol');
hs=findobj('Tag','mfload_scan');
hn=findobj('Tag','mfload_norm');
hp=findobj('Tag','mfload_pol');
%hf1=findobj('Tag','mfload_F1');
%hf2=findobj('Tag','mfload_F2');

if ~isempty(scan) scano = scan; else scano = []; end
if isempty(xcol)   xcol=get(hx,'Value'); end
if isempty(ycol)   ycol=get(hy,'Value'); end
if isempty(ecol)   ecol=get(he,'Value'); end
if isempty(mcol)   mcol=get(hm,'Value'); end
if isempty(scan)   scan=get(hs,'String'); end
if isempty(normi)  normi=get(hn,'String'); end
if isempty(xcol) | (xcol > ncol+1), xcol = ncol+1; end
if isempty(ycol) | (ycol > ncol), ycol = 1; end
if isempty(ycol) | (ecol > ncol+2), ecol = ncol+1; end
mcol = min(mcol,ncol+1);

set(hx,'String',[collist '|Point Number'],'Value',xcol(1));
set(hy,'String',collist,'Value',ycol(1));
set(he,'String',[collist '|sqrt(y)|none'],'Value',ecol(1));
set(hm,'String',[collist '|none'],'Value',mcol(1));
set(hs,'String',scan);
set(hn,'String',normi);
if ~isempty(Fstate)
  if (get(hp,'Value') > max(size(Fstate)))
    set(hp,'Value',1)
  end
  set(hp,'String',plist);
end

disp('REM: the ''scan'' field can be used to enter options such as ''1,setpar''.');

set(findobj('Tag','mf_coldlg'),'Visible','on');          % Switch window on
delete(gca);

%============ Now wait until ok or cancel pressed =========================
scano = scan;
drawnow;
hok = findobj('Tag','mfload_ok');
hcan= findobj('Tag','mfload_cancel');
my_gco = hok;

if isempty(file_index)
  % should only wait if not in a serie or column types changed

  if(0)
    waitforbuttonpress;

    while (~isempty(gco) & (gco ~= hok) & (gco ~= hcan))
        drawnow

    end
    my_gco = gco;

  end


  if(1)
    while (isempty(gco))
        %pause(1)
        drawnow
    end


%     while( ~isempty(gco) & (gco ~= hok) & (gco ~= hcan))
    while((hok ~= gco) & (hcan ~=gco ))
        %pause(1)
        drawnow
    end

    if(isempty(gco))
        axel
    end

    my_gco = gco;
    delete(gco);

    uicontrol(hmf_load,...
       'Style','PushButton',...
       'String','Ok',...
       'Tag','mfload_ok',...
       'Position',[32 4 50 20]);

  end
end

%========= Extract relevant data columns ==================================
if (~isempty(my_gco) & (my_gco == hok))

%---------- Retrieve column numbers and labels from figure ---------------
   collist=get(findobj('Tag','mfload_xcol'),'String');      % (in nice mat string format)
   xcol=get(hx,'Value');          % menu positions
   ycol=get(hy,'Value');
   ecol=get(he,'Value');
   mcol=get(hm,'Value');
   scan=get(hs,'String');
   normi=get(hn,'String');
   flip=get(hp,'Value');
   xlab=collist(xcol,:);
   ylab=collist(ycol,:);


  scan = num2str(scan);
  scano = num2str(scano);

   t1 = sscanf(scano,'%i');
   t2 = sscanf(scan,'%i');
   if isempty(t2) | (~isempty(t1) & ( t1 == t2 ))
    set(findobj('Tag','mf_coldlg'),'Visible','off');
   else
    keep = 1;
   end
   xcolo =xcol;
   ycolo=ycol;
   mcolo=mcol;
   ecolo=ecol;
   scano=scan;
   normo=normi;
   flipo=flip;

else
   xcol=[]; ycol=[]; ecol=[]; mcol=[]; xlab=''; ylab=''; scano=''; normo=''; keep = 0; flip=[0,0];
   delete(findobj('Tag','mf_coldlg'));
end


