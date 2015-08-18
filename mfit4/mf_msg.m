function mf_msg(msg)
%
% MFIT  function  mf_msg(msg)
%     Display message in message window
%     MZ 29.11.94
%
set(findobj('Tag','mf_MessageWindow'),'String',msg);
fprintf(1, 'MFIT: %s\n', msg);

