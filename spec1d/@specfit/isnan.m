function [ output_args ] = isnan( input_args )
%ISNAN Summary of this function goes here
%   Detailed explanation goes here
output_args = 0;

if any(isnan(input_args.pvals))
output_args = 1;
end

