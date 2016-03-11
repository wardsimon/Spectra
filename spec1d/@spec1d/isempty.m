function [ output_args ] = isempty(varargin)
%ISEMP Summary of this function goes here
%   Detailed explanation goes here
j = 1;
for i = 1:length(varargin)
    objs = varargin{i};
    for k = 1:length(objs);
        if length(objs(k).x) == 0
            output_args(j) = 1;
        else
            output_args(j) = 0;
        end
        j = j+1;
    end
end
end

