function ml_zoom
%
% function ml_zoom
%
% MLAUE function to control zoom in on and extracting data
% ARW 15.10.06
%
% Last modified: ARW 16.10.06

%===== Find the window ============================================
hml_ctrl = findobj(0,'Tag','ml_ControlWindow'); 
if isempty(hml_ctrl) hml_ctrl = 0; end

%===== Find the status of the buttons =============================
hml_zoom=findobj(hml_ctrl,'Tag','ml_ZoomRadio');
butstat=[get(hml_zoom,'Value')];
if butstat(1,1)==1
     zoom on;
else
     zoom off;
end

