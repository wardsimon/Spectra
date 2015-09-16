function ml_extplot
%
% function ml_extplot
%
% MLAUE function to plot the extracted data
% ARW 05.08.07
%
% Last modified:

%===== Find the windows ============================================
hml_ctrl = findobj(0,'Tag','ml_ControlWindow'); 
if isempty(hml_ctrl) hml_ctrl = 0; end
hml_imwin = findobj(hml_ctrl,'Tag','ml_ImageWindow');
hml_exwin = findobj(0,'Tag','ml_ExtDataWindow');
if isempty(hml_exwin) hml_exwin = 0; end
hml_exdat=findobj('Tag','ml_ExtractedImage');
hml_bkgrad=findobj('tag','ml_BackGroundRadio');

%----- Get the extracted data, stored as UserData in the Extracted Data Window ---
exdat=get(hml_exdat,'UserData');
if isempty(exdat), exdat=get(hml_exwin,'UserData');end
%----- Get the background borders, stored as UserData in the Select Background Button -------
borders=get(hml_bkgrad,'UserData');

%----- Invert the limits defined in the main MLAUE window ---------------------
imclim=65535-fliplr(get(hml_imwin,'Clim'));
%----- Plot data -----
axes(hml_exdat);
imagesc(exdat,imclim);
%----- Use the inverse colormap as defined in the main MLAUE window -----------
cmap=get(findobj(hml_ctrl,'Tag','ml_ImageMenu'),'UserData');
colormap(hml_exdat,eval(['flipud(',cmap,')']));
%----- Plot a Background Box if it is defined --------------------------------
if ~isempty(borders)
    hold on;
    nbox=length(borders)/5;
    for i=1:nbox
        for j=1:4
            index=5*(i-1)+j;
            plot([borders(index,1),borders(index+1,1)],[borders(index,2),borders(index+1,2)]);
        end
    end
    hold off;
end
set(hml_exdat,'Tag','ml_ExtractedImage','UserData',exdat,'Visible','On');
