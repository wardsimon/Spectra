function ml_modplot
%
% function ml_modplot
%
% MLAUE function to modify the limits of the image
% ARW 15.10.06
%
% Last modified:ARW 13.8.07

%===== Find the window ============================================
hml_ctrl = findobj(0,'Tag','ml_ControlWindow'); 
if isempty(hml_ctrl) hml_ctrl = 0; end

%===== Find the image =============================================
hml_imwin = findobj(hml_ctrl,'Tag','ml_ImageWindow');
axes(hml_imwin);

%===== Find the min/max parameters ================================
Imin=sscanf((get(findobj('tag','ml_MinIntens'),'String')),'%g');
Imax=sscanf((get(findobj('tag','ml_MaxIntens'),'String')),'%g');

%===== Check that min is not > max ================================
if Imin >= Imax
    disp('ERROR: The minimum must be greater than or equal to the maximum');
    return
end

%===== Determine whether the image fix is set =====================
hifrb=findobj('tag','ml_IntFix');
if get(hifrb,'Value')==1
    set(hifrb,'UserData',[Imin,Imax]);
else
    set(hifrb,'UserData',[]);
end
    
%===== Set the min/max parameters =================================
set(gca,'Clim',[Imin,Imax]);