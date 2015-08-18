function ml_msg(txt)
%
% function ml_msg(txt)
%
% MLAUE function to send text to the extracted data window
% ARW 03.08.07
%
% Last modified:

%===== Find the extracted data window ============================================
hml_exwin = findobj(0,'Tag','ml_ExtDataWindow'); 
if isempty(hml_exwin) hml_exwin = 0; end
%===== Find the message window ============================================
hml_msg = findobj(hml_exwin,'Tag','ml_msg');

set(hml_msg,'string',txt);
