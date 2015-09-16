function ml_load(datadir,datafile);
%
% function ml_load(datadir,datafile);
%
% MLAUE function to load and display a TIFF image
% ARW 14.10.06
%
% Last modified: ARW 13.8.07

lims=[];

%===== Find the window ============================================
hml_ctrl = findobj(0,'Tag','ml_ControlWindow'); 
if isempty(hml_ctrl) hml_ctrl = 0; end

%===== Find the data file/directory ===============================
if nargin ~=2,
    datadir=get(findobj('Tag','ml_DataDir'),'string');           % data directory
    if isempty(datadir), datadir = pwd; end
    datafile=get(findobj('Tag','ml_DataFile'),'string');         % data file
    if isempty(datafile), [datafile, datadir]=uigetfile('*.tif'); end
    lims=get(findobj('tag','ml_IntFix'),'UserData');
end

%===== Get the data ===============================================
data=imread([datadir,datafile],'tif');
if isempty(lims), lims=[min(min(data)) max(max(data))]; end

%===== Display the data ===========================================
hml_imwin = findobj(hml_ctrl,'Tag','ml_ImageWindow');
hml_imwin_flags=get(hml_imwin);
axes(hml_imwin);
imagesc(data,lims);
colormap(gray);
set(gca,'Tag','ml_ImageWindow');
set(gca,'UserData',data);

%===== Update the text ============================================
set(findobj('Tag','ml_DataDir'),'String',datadir);
set(findobj('Tag','ml_DataFile'),'String',datafile);
set(findobj('Tag','ml_MinIntens'),'String',lims(1));
set(findobj('Tag','ml_MaxIntens'),'String',lims(2));

%===== Set the zoom/extract buttons to 'off' ======================
set(findobj(hml_ctrl,'Tag','ml_ZoomRadio'),'Value',0);
zoom off;
set(findobj(hml_ctrl,'Tag','ml_ExtractRadio'),'Value',0);
