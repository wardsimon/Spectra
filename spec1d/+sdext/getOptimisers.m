function [ output_args ] = getOptimisers()
%GETOPTIMISERS
% Get a list of available optimisers for fitting.

d = fullfile(sdext.getpref('libroot').val,'Spectra','spec1d','Libraries','Optimizers','*.m');
files = dir(d);
output_args = strrep({files.name},'.m','');
end

