%% @spec1d/plot
% This is the reference documentation for the function @spec1d/plot
%
% This function creates a plot of a given spectra or vecor of spec1d objects. 
%%

%% Syntax
%
%    varagout = plot(s,varargin)

%% Inputs
%
% * _s_ - Single or vector of spec1d objects
% 

%% Optional Inputs
%
% There are a few built in common options. Which can be used by passing
% parameter/value pairs. These are:
%
% The option *semilogy*
% * *0*: The y-axis is not in log scale. (Default)
% * *1*: The y-axis is not in log scale.
%
% The option *semilogx*
% * *0*: The x-axis is not in log scale. (Default)
% * *1*: The x-axis is not in log scale.
%
% The option *loglog*
% * *0*: The x-axis and y-axis is not in log scale. (Default)
% * *1*: The x-axis and y-axis  is not in log scale.
%
% The option *tLength* sets the length of the error bar tees
% * *NaN*: tee length is calculated as 0.075 * the minimum distance between points.
% * *x>0*: tee length is calculated as *x* * the minimum distance between points. 
% * *x<0*: tee length is supplied by the value *x*.
%

%% Optional Outputs
%
% * _s_out_ - Spec1d object of the same size as the input array

%% Example
% This is an example on using @spec1d/plot
%
% <html><h3>Example 1</h3></html>
%

s = spec1d(1:0.2:5,sin(linspace(0,2*pi,21)),0.1);
s_1 = plot(s);

figure
plot(s,s_1)
legend({'s','plot(s)'})

%% See Also
% <html><a href="{{ site.url }}/@spec1d.plus/index.html">Plus</a>, <a href="{{ site.url }}/@spec1d.sum/index.html">Sum</a></html>