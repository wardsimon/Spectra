%% @spec1d/abs 
% This is the reference documentation for the function @spec1d/abs
%
% This function gives a new spectra where the absolute y values of the
% imput spectra _s_ are used.
%%

%% Syntax 
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
% 
% <html><h3>Example 1</h3></html>
%

s = spec1d(1:0.2:5,sin(linspace(0,2*pi,21)),0.1);
s_1 = abs(s);

figure
plot(s,s_1)
legend({'s','abs(s)'},'Location','SouthWest')