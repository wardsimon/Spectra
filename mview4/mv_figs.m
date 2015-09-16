function [hmv_ctrl, hmv_data]=mv_figs
%
% Mview function [hmv_ctrl, hmv_data]=mv_figs
% 	Returns handles of control, data, and figure windows;
% 	zero if they don't exist
% 	MZ 29.11.94

CtrlName='MVIEW: Control';
DataName='MVIEW: Data';
hmv_ctrl=0;
hmv_data=0;

objs = get(0,'Children');
for i = 1:length(objs),
   if strcmp(get(objs(i),'type'),'figure')
      Name=get(objs(i),'Name');
      if strcmp(Name,CtrlName) hmv_ctrl=objs(i); end;
      if strcmp(Name,DataName) hmv_data=objs(i); end;
   end 
end


