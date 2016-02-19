%% @spec1d/set
% This is the reference documentation for the function @spec1d/set
%
% This function is an implementation of the get/set routine from the handle
% class. Using _set_ it is possible to set the field/s in a spec1d object.
%%

%% Syntax
%
%    s_out = set(s,'m',m,...)
%

%% Inputs
%
% * _s_ - Singlespec1d object.
% * '_m_' - Fieldname to be set.
% * _m_ - Vector corresponding to field '_m_'.
%

%% Outputs
%
% * _s_out_ - Spec1d object with changed field/s

%% Example
% This is an example on using @spec1d/set
%
% <html><h3>Example 1</h3></html>
%

s = spec1d(1:0.2:5,sin(linspace(0,2*pi,21)),0.1);
s_1 = set(s,'e',0.1*rand(21,1));

fprintf('e (s)\te (s_1)\n')
fprintf('%0.03f\t%0.03f\n',[get(s,'e') get(s_1,'e')]')

%% Note
% * The vector being set must have an equal number of points as the current fields in _s_.  
% * Multiple fields in _s_ can be changed at once by adding arguments.

%% See Also
% <html><a href="{{ site.url }}/@spec1d.get/index.html">Get</a></html>