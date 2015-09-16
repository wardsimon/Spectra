function mf_fitrt
%
% MFIT function mf_fitrt
%     Load new fitting routine
%     MZ 29.11.94
%
[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

%---------- Get fit function details ---------------
h=get(hmf_ctrl,'userdata');   % handles

%---------- Choose fit routine -------------------
err=mf_file('Select fit routine:',h(10),h(11),'*.m');
if err==0
  fitdir=get(h(11),'String');
  fitfile=get(h(10),'String');
  fitfile=lower(fitfile(1:findstr(fitfile,'.')-1));
  set(h(10),'String',fitfile);

  %----------- Update func names ----------------------------
  set(h(9),'String',fitfile);       % Update func name
  addpath(fitdir);
end
