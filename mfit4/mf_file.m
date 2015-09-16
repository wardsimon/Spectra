function err=mf_file(titl, hfile, hdir,filter)
%
% MFIT function  err=mf_file(titl, hfile, hdir)
%     Dialog to choose filename
%     MZ 29.11.94

if nargin==3
   filter='*';
end

%---------- Get dir name --------------
newdir=get(hdir,'String');

%---- Dialogue to choose file -------------
[file, fdir]=uigetfiles([newdir filter],titl,'MultiSelect','off');

if file~=0
  set(hfile,'String',file);
  set(hdir,'String',fdir);
  err=0;
else
  err=1;
end;


