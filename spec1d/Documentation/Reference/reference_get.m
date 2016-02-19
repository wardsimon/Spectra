%% @spec1d/get
% This is the reference documentation for the function @spec1d/get
%
% This function is an implementation of the get/set routine from the handle
% class. Using _get_ it is possible to get the fields in a spec1d object.
%%

%% Syntax
%
%    m = get(s,'m');
%    [m, n] = get(s,'m','n');
%

%% Inputs
%
% * _s_ - Single or vector of spec1d objects
%

%% Outputs
%
% * _m_ - Column vecor of field '_m_'.

%% Examples
% These are examples on using @spec1d/get
%
% <html><h3>Example 1</h3></html>
%
% Get the x-values of a spectra
%

s = spec1d(1:0.2:5,sin(linspace(0,2*pi,21)),0.1);
x = get(s,'x');
fprintf('x\n')
fprintf('%0.03f\n',x')

%% 
% <html><h3>Example 2</h3></html>
%
% Get the x and y values of a spectra
%

[x,y ] = get(s,'x','y');
fprintf('x\ty\n')
fprintf('%0.03f\t%0.03f\n',[x y]')

%% See Also
% <html><a href="{{ site.url }}/@spec1d.extract/index.html">Extract</a>, <a href="{{ site.url }}/@spec1d.set/index.html">Set</a></html>