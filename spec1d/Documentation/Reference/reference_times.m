%% @spec1d/times 
% This is a the reference documentation for the function @spec1d/times
% This function gives a new spectra with a value of spectra _s1_ multiplied by a value or spectra _s2_.
%%

%% Usage 
% 
%  # s_out = times(s1,s2)
%  # s_out = s1 * s2
%  # s_out = s1 * a
%  # s_out = a * s1

%% Inputs
% 
% * _s1_ - spec1d object
% * _s2(a)_ - spec1d object or value
% 

%% Outputs
% 
% * _s_out_ - spec1d object

%% Notes
% 
% * If spectra are not the same length, _s2_ is interpolated to the x
% values of _s1_.
% 

%% Examples
% These are some examples on using @spec1d/times

%%% Example 1
% Multiplication of spectra and value; 

s = spec1d(1:5,1:5,0.1);
s_1 = s * 5;
s_1 = 5 * s;

%%% Example 2
% Multiplication two spectra (equal points)

s2 = spec1d(1:5,2:6,0.15);
s_2 = s * s2;

%%% Example 3
% Addition two spectra (unequal points). This using interpolation on _s3_

s3 = spec1d(1:5 + 0.1,2:6,0.15);
s_3 = s * s3;
