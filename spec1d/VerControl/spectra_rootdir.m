function rootdir = spectra_rootdir()
% gives the path to the Spectra code
%
% rootdir = SPECTRA_ROOTDIR()
%
% See also SPECTRA.
%

rootdir = mfilename('fullpath');
idx     = strfind(rootdir,filesep);
rootdir = rootdir(1:idx(end-2));

end