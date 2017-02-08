function rootdir = spectra_rootdir()
% gives the path to the Spectra code
%
% rootdir = SPECTRA_ROOTDIR()
%
% See also SPECTRA.
%
% This has been modified from the spinW help file. Thanks to Sandor Toth!

rootdir = mfilename('fullpath');
idx     = strfind(rootdir,filesep);
rootdir = rootdir(1:idx(end-2));

end