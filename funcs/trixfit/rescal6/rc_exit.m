% Rescal function rc_exit
% Clean exit from Reciprocal
% D.A.T adapted from M.Z. 16.10.95
% ========= Close the windows ===========================
disp('Closing windows');

delete(findobj('tag','Rescal: Parameters'));
delete(findobj('name','Rescal: Projections'));
delete(findobj('tag','Rescal: Instrument'));
delete(findobj('tag','Rescal: Reciprocal'));
delete(findobj('tag','Rescal: Spectrometer'));

% ===== Done ============================================

t=fix(clock);
disp(sprintf('Rescal ended. %d.%d.%d.   %d:%d:%d\n',t(3:-1:1),t(4:6)));

