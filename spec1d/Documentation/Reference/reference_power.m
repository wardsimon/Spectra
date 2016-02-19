%% @spec1d/power
% This is the reference documentation for the function @spec1d/power
%
% This function raises the y-values of a single or vecor of spec1d
% objects to a power specified.
%%

%% Syntax
%
%    s_out = power(s,m)
%    s_out = s^m
%    s_out = s.^m
%

%% Inputs
%
% * _s_ - Single or vector of spec1d objects.
% * _m_ - Single value corresponding to the power to raise spectra _s_ by.
%

%% Outputs
%
% * _s_out_ - Spec1d object of the same size as the input array

%% Example
% This is an example on using @spec1d/power
%
% <html><h3>Example 1</h3></html>
%

s = spec1d(1:0.2:5,linspace(1,2,21),0.1);
s_1 = s.^2;

figure
plot(s,s_1)
legend({'s','power(s,2)'},'Location','NorthWest')

%% See Also
% <html><a href="{{ site.url }}/@spec1d.sqrt/index.html">Sqrt</a></html>