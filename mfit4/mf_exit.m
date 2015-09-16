function  mf_exit(cmd)
% MFIT function mf_exit
% Clean exit from MFIT
%
% M. Zinkin 8.3.95

%=====Delete MFIT path from matlab path =============================

%disp('Removing Mfit directories from path');

[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

% autosave for config ============================================
h=findobj('tag','mf_AutoSave');
if ~isempty(h) & ~isempty(get(h,'string'))
  t=fix(clock);
  disp('Closing windows and removing Mfit objects');
  disp(sprintf('Mfit end. %d.%d.%d   %d:%d:%d\n',t(3:-1:1),t(4:6)));
  asv = str2num(get(h,'string'));
  if ~isempty(asv) & (asv == 1)
    mf_save('ini');
    set(h,'string',[])
  end
end

%===== Close the windows ============================================
hmat = get(0,'children');
for i = 1:size(hmat)
  h = hmat(i);
  tag = eval('get(h,''tag'')','[]');
  if (h~=hmf_ctrl) & (~isempty(strmatch('mf',tag)) | ~isempty(strmatch('hmf',tag)))
    delete(h);
  end
end

if ~isempty(dir('mf_mftmp.m'))
  disp('Note : Mfit multifunction temporary file mf_mftmp.m removed.')
  delete('mf_mftmp.m');
end

delete(hmf_ctrl);

