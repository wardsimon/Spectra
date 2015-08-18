function ml_reset

% function ml_reset
% MLAUE function to clear all background, mask information 
% and replot the original extracted data
%
% ARW 11.08.07

%===== Find all the relevant objects =====
hml_exwin=findobj('Tag','ml_ExtDataWindow');    % Contains the original data
hml_exdat=findobj('Tag','ml_ExtractedImage');   % Contains the (masked) data
hml_bkgtit=findobj('tag','ml_BkgTit');          % Contains the x,y coordinates for the background
hml_bkgrad=findobj('tag','ml_BackGroundRadio'); % Contains bounds for the background 

%===== Get the original data =====
data=get(hml_exwin,'UserData');

%===== Reset the UserData fields =====
set(hml_exdat,'UserData',data);
set(hml_bkgtit,'UserData',[]);
set(hml_bkgrad,'UserData',[]);

%----- Calculate Integrated Intensities ---------------------------------------
ml_updint;
%----- Plot data -----
ml_extplot;
