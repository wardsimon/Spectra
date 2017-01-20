function [ output_args ] = getCriteria()
%GETCRITERIA
% Get a list of available criteria for fitting.

d = fullfile(sdext.getpref('libroot').val,'Spectra','spec1d','Libraries','Criteria','*.m');
files = dir(d);
output_args = strrep({files.name},'.m','');
end


