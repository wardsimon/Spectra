function  mf_exit
% Mview function mv_exit
% Clean exit from Mview
%
% EF 03.98

%===== Close the windows ============================================
disp('Closing windows and removing Mview objects');
hmat = get(0,'children');
for i = 1:size(hmat)
  h = hmat(i);
  try
    tag = lower(get(h,'tag'));
  catch
    tag='';
  end
  if ~isempty(tag) & ~isempty(h)
  if ~isempty(strmatch('mv',tag)) | ~isempty(strmatch('hmv',tag)) | ~isempty(strmatch('mf_coldlg',tag)) | ~isempty(strmatch('mf_cb',tag))
    delete(h);
  end
  end
end

%===== Done =========================================================
t=fix(clock);
disp(sprintf('Mview end. %d.%d.%d   %d:%d:%d\n',t(3:-1:1),t(4:6)));

