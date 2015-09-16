function p = iInstall
% iInstall: installs the iLib (ILL Library) directories
%
% Part of: iLoad utilities (ILL library)
% Author:  E. Farhi <farhi@ill.fr>. Feb 4th, 2003.

local_path = fileparts(which('iInstall'));
to_add = {'.', 'iFiles', 'iStrings', 'iLoad', 'iGrep'};
for index=1:length(to_add)
  addpath(fullfile(local_path, '..', to_add{index}));
end
disp('iLib: Installed the ILL library directories.');
