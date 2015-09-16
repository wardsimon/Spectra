function [hmf_ctrl, hmf_data, hmf_pars]=mf_figs
%
% MFIT function [hmf_ctrl, hmf_data, hmf_pars]=mf_figs
% 		Returns handles of control, data, and figure windows;
% 		zero if they don't exist
% 		MZ 29.11.94

CtrlName = 'mf_ControlWindow';
DataName = 'mf_DataWindow';
ParsName = 'mf_ParWindow';

hmf_ctrl = findobj(0,'Tag',CtrlName); 
if isempty(hmf_ctrl) hmf_ctrl = 0; end

hmf_data = findobj(0,'Tag',DataName); 
if isempty(hmf_data) hmf_data = 0; end

hmf_pars = findobj(0,'Tag',ParsName); 
if isempty(hmf_pars) hmf_pars = 0; end
