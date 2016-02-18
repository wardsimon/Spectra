%% @spec1d/cumsum 
% This is the reference documentation for the function @spec1d/cumsum
%
% This function gives the cumulative summation of points in each spectra _s_ were _s_ can be single or an array.
%%

%% Syntax 
% 
%    s_out = cumsum(s)

%% Inputs
% 
% * _s_ - Single or vector of spec1d objects 
% 

%% Outputs
% 
% * _s_out_ - Spec1d object of the same size as the input array

%% Example
% This is an example on using @spec1d/cumsum
% 
% <html><h3>Example 1</h3></html>
%

s = spec1d(1:0.2:5,sin(linspace(0,2*pi,21)),0.1);
s_1 = cumsum(s);

figure
plot(s,s_1)
legend({'s','cumsum(s)'})

%% See Also
% <html><a href="{{ site.url }}/@spec1d.plus/index.html">Plus</a>, <a href="{{ site.url }}/@spec1d.sum/index.html">Sum</a></html>