function ml_exwin
%
% function ml_exwin
%
% MLAUE function to open a window to manipulate extracted data
% from the image
%
% ARW 16.10.06
% Last modified:  13.08.07

%===== Find the windows ============================================
hml_ctrl = findobj(0,'Tag','ml_ControlWindow'); 
if isempty(hml_ctrl) hml_ctrl = 0; end
hml_imwin = findobj(hml_ctrl,'Tag','ml_ImageWindow');
hml_exwin = findobj(0,'Tag','ml_ExtDataWindow');
if isempty(hml_exwin) hml_exwin = 0; end
%===== Switch the ZOOM off ========================================
set(findobj('tag','ml_ZoomRadio'),'Value',0);ml_zoom;

%===== Extract the data and convert to double ======================
data=get(hml_imwin,'UserData');
datlims=[get(hml_imwin,'Ylim') get(hml_imwin,'Xlim')];
exdat=65535-double(data(datlims(1):datlims(2),datlims(3):datlims(4)));

%===== Set the window sizes ========================================
wpos=[ 200 20 800 350 ];
expos=[0.55 0.05 0.425 0.90];
Vpos=[25 5];
Hpos=[10 115];

%===== Create the window ===========================================
if hml_exwin ~= 0
    figure(hml_exwin);
    hml_exdat=findobj('Tag','ml_ExtractedImage');
    axes(hml_exdat);
else
    hml_exwin = figure('Position',wpos,...
            'tag','ml_ExtDataWindow',...
            'Name','MLAUE: Extracted Data',...
            'DefaultUicontrolBackgroundColor',[1 1 1],...
            'DefaultUicontrolForegroundColor',[0 0 0],...
            'DefaultUicontrolHorizontalAlignment','left',...
            'ToolBar','figure',...
            'MenuBar','none',...
            'visible','on',...
            'HandleVisibility','on',...
            'NextPlot','Add',...
            'Interruptible','on');
%----- Message box -----
    hml_msg=uicontrol(hml_exwin,...
            'style','text',...
            'string','Extracted Data Manipulation Window',...
            'Tag','ml_msg',...
            'position',[Hpos(1) wpos(4)-sum(Vpos) 3*Hpos(2) Vpos(1)]);    

%----- Reset image -----
    uicontrol(hml_exwin,...
        'style','pushbutton',...
        'String','Clear and reset',...
        'callback','ml_reset',...
        'Position',[Hpos(1) wpos(4)-2*sum(Vpos) Hpos(2) Vpos(1)]);
%----- Radiobutton to mask points -----
    hml_mskrad=uicontrol(hml_exwin,...
            'tag','ml_MaskRadio',...
            'style','radiobutton',...
            'String','Mask Points',...
            'callback','ml_mskpt',...
            'Value',0,...
            'UserData',0,...
            'BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'), ...
            'ForegroundColor','k',...
            'Position',[Hpos(1) wpos(4)-4*sum(Vpos) Hpos(2) Vpos(1)]);
%----- Radiobutton to choose background -----
    hml_bkgrad=uicontrol(hml_exwin,...
            'tag','ml_BackGroundRadio',...
            'style','radiobutton',...
            'String','Select Background',...
            'callback','ml_bkgd',...
            'Value',0,...
            'UserData',[],...
            'BackgroundColor',get(0,'DefaultUicontrolBackgroundColor'), ...
            'ForegroundColor','k',...
            'Position',[Hpos(1) wpos(4)-5*sum(Vpos) Hpos(2) Vpos(1)]);

%===== Integrated data windows =====
%----- Integrated Counts -----------------
    uicontrol(hml_exwin,...
        'style','text',...
        'String','Integrated Counts',...
        'Position',[Hpos(1) wpos(4)-6*sum(Vpos) Hpos(2) Vpos(1)]);
    hml_TotIntBox=uicontrol(hml_exwin,...
        'style','text',...
        'Tag','ml_TotIntBox',...
        'UserData',0,...
        'Position',[Hpos(1)+sum(Hpos) wpos(4)-6*sum(Vpos) Hpos(2) Vpos(1)]);
    uicontrol(hml_exwin,...
        'style','text',...
        'String','+/-',...
        'Position',[Hpos(1)+2*sum(Hpos) wpos(4)-6*sum(Vpos) Hpos(1) Vpos(1)]);
    hml_eTotIntBox=uicontrol(hml_exwin,...
        'style','text',...
        'Tag','ml_eTotIntBox',...
        'UserData',0,...
        'Position',[3*Hpos(1)+2*sum(Hpos) wpos(4)-6*sum(Vpos) Hpos(2) Vpos(1)]);
%----- Background -----------------
    uicontrol(hml_exwin,...
        'style','text',...
        'String','Background per point',...
        'tag','ml_BkgTit',...
        'UserData',[],...
        'Position',[Hpos(1) wpos(4)-7*sum(Vpos) Hpos(2) Vpos(1)]);
    hml_BkgBox=uicontrol(hml_exwin,...
        'style','text',...
        'Tag','ml_BkgBox',...
        'UserData',0,...
        'Position',[Hpos(1)+sum(Hpos) wpos(4)-7*sum(Vpos) Hpos(2) Vpos(1)]);
    uicontrol(hml_exwin,...
        'style','text',...
        'String','+/-',...
        'Position',[Hpos(1)+2*sum(Hpos) wpos(4)-7*sum(Vpos) Hpos(1) Vpos(1)]);
    hml_eBkgBox=uicontrol(hml_exwin,...
        'style','text',...
        'Tag','ml_eBkgBox',...
        'UserData',0,...
        'Position',[3*Hpos(1)+2*sum(Hpos) wpos(4)-7*sum(Vpos) Hpos(2) Vpos(1)]);
%----- Integrated Counts - Background -----------------
    uicontrol(hml_exwin,...
        'style','text',...
        'String','Int - Bkg',...
        'Position',[Hpos(1) wpos(4)-8*sum(Vpos) Hpos(2) Vpos(1)]);
    hml_IminusBBox=uicontrol(hml_exwin,...
        'style','text',...
        'Tag','ml_IminusBBox',...
        'UserData',0,...
        'Position',[Hpos(1)+sum(Hpos) wpos(4)-8*sum(Vpos) Hpos(2) Vpos(1)]);
    uicontrol(hml_exwin,...
        'style','text',...
        'String','+/-',...
        'Position',[Hpos(1)+2*sum(Hpos) wpos(4)-8*sum(Vpos) Hpos(1) Vpos(1)]);
    hml_eIminusBBox=uicontrol(hml_exwin,...
        'style','text',...
        'Tag','ml_eIminusBBox',...
        'UserData',0,...
        'Position',[3*Hpos(1)+2*sum(Hpos) wpos(4)-8*sum(Vpos) Hpos(2) Vpos(1)]);
%----- Create image axes for extracted data  --------
    hml_exdat=axes('Position',expos,...
                   'Tag','ml_ExtractedImage',...
                   'UserData',exdat);

end

%----- Invert the limits defined in the main MLAUE window ---------------------
imclim=65535-fliplr(get(hml_imwin,'Clim'));
%----- Store the original data in case of error -------------------------------
set(hml_exwin,'UserData',exdat);
%----- Calculate Integrated Intensities ---------------------------------------
ml_updint;
%----- Plot data -----
ml_reset;