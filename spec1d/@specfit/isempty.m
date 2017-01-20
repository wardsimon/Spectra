function [ output_args ] = isempty( obj )
%ISEMPTY Summary of this function goes here
%   Detailed explanation goes here
output_args = zeros(1,length(obj));
for i = 1:length(obj)
    if isempty(obj(i).pvals)
        output_args(i) = 1;
    end
end

