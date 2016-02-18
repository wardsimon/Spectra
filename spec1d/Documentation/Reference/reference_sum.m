%% @spec1d/sum 
% This is the reference documentation for the function @spec1d/sum
%
% This function gives the summation of points in each spectra _s_ where _s_
% can be single or an array. The error _ee_ can also be optionally
% returned.
%%

%% Syntax 
% 
%    yy = sum(s)
%    [yy, ee] = sum(s)

%% Inputs
% 
% * _s_ - Single or vector of spec1d objects 
% 

%% Outputs
% 
% * _yy_ - Vector of length _s_
% * _ee_ - Vector of length _s_ (optional)

%% Examples
% This is an example on using @spec1d/sum
%
% <html><h3>Example 1</h3></html>
%
% Return the sum of each sprecta

s(1) = spec1d(1:5,1:5,0.1);
s(2) = spec1d(1:5,2:6,0.075);

yy = sum(s);
disp(yy)

%%
% 
% <html><h3>Example 2</h3></html>
%
% Return the sum of each sprecta with error
%

[yy, ee] = sum(s);
fprintf('s(1)\ts(2)\n')
fprintf('%0.03f\t%0.03f\n',yy)
fprintf('%0.03f\t%0.03f\n',ee)

%% See Also
% <html><a href="{{ site.url }}/@spec1d.plus/index.html">Plus</a>, <a href="{{ site.url }}/@spec1d.cumsum/index.html">Cumsum</a></html>