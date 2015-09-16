function mv_msg(msg)
%
% Mview  function  mv_msg(msg)
% Display message in message window
% DFM 5.6.96
%
h=findobj('Tag','hmv_message');
if ~isempty(h); set(h,'String',msg); end


