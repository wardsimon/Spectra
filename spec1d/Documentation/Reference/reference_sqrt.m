%% @spec1d/sqrt
% This is the reference documentation for the function @spec1d/sqrt
%
% This function takes the square root of a single or vecor of spec1d
% objects.
%%

%% Syntax
%
%    s_out = sqrt(s)

%% Inputs
%
% * _s_ - Single or vector of spec1d objects
%

%% Outputs
%
% * _s_out_ - Spec1d object of the same size as the input array

%% Example
% This is an example on using @spec1d/sqrt
%
% <html><h3>Example 1</h3></html>
%

s(1) = spec1d(1:5,1:5,0.1);
s(2) = spec1d(1:0.2:5,sin(linspace(0,2*pi,21)),0.1);
s_1 = sqrt(s(1));

figure
plot(s(1),s_1)
legend({'s','sqrt(s)'},'Location','NorthWest')


%%
%
% <html><h3>Example 2</h3></html>
%
% Immaginary numbers are not supported
%

s_2 = sqrt(s(2));

%% See Also
% <html><a href="{{ site.url }}/@spec1d.power/index.html">Power</a></html>