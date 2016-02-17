%% @spec1d/abs 
% This is a the reference documentation for the function @spec1d/abs
%
% This function gives a new spectra where the absolute y values of the
% imput spectra _s_ are used.
%%

%% Usage 
% 
%    s_out = abs(s)

%% Inputs
% 
% * _s_ - single or vector of spec1d objects 
% 

%% Outputs
% 
% * _s_out_ - spec1d object of the same size as the input array

%% Example
% This is an example on using @spec1d/abs

%% Example 1

s = spec1d(-5:5,-5:5,0.1);
s_1 = abs(s);
