function  ml_exit(cmd)
% MLAUE function ml_exit
% Clean exit from MLAUE
%
% ARW 14.10.06
% Based on MF_EXIT by M. Zinkin 8.3.95
%
% Last modified:

%===== Find the windows ============================================
hml_ctrl = findobj(0,'Tag','ml_ControlWindow'); 
if isempty(hml_ctrl) hml_ctrl = 0; end

%===== Close the windows ============================================
hmat = get(0,'children');
for i = 1:size(hmat)
  h = hmat(i);
  tag = eval('get(h,''tag'')','[]');
  if (h~=hml_ctrl) & (~isempty(strmatch('ml',tag)) | ~isempty(strmatch('hml',tag)))
    delete(h);
  end
end

delete(hml_ctrl);

