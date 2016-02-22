%% @spec1d/spec1d
% This is the reference documentation for the function @spec1d/spec1d
%
% This function creates a spec1d object.
%%

%% Syntax
%
%    s_out = spec1d(x,y,e,varargin);
%    s_out = spec1d(struct);
%    s_out = spec1d(s);
%

%% Inputs
%
% Create a spec1d object with data points.
% * _x_ - Single or vector of x-points for the spec1d object
% * _y_ - Single or vector of y-points for the spec1d object
% * _e_ - Single or vector of e-points for the spec1d object
%
% Create a spec1d object with a structure
% * _struct_ - Structure with fields x, y, e.
%
% Reform or validate a spec1d object.
% * _s_ - A spec1d object.
%

%% Optional Inputs
% * _x_label_ - String, automatic x-label for plotting.
% * _y_label_ - String, automatic y-label for plotting.
% * _datafile_ - String, automatic title for plotting or data identifier.
% * _yfit_ - Single or vector of y-points corresponding to a fit. 

%% Outputs
%
% * _s_out_ - Spec1d object

%% Notes
% * Error can be given as a number and replicated for all y-points.  
% * If no arguments are given an empty spec1d is returned.
% * If not all fields are present, they are padded with zeros. For example
% spec1d(0:5) will retrun a spec1d object with x = 0:5 and with zeros in y
% and e.
% * Pre 2013 spec1d objects are converted into the new data structure
% automatically.

%% Examples
% These are examples on using @spec1d/spec1d
%
% <html><h3>Example 1</h3></html>
%

x = 1:0.2:5;
y = sin(linspace(0,2*pi,21));
e = 0.1;

s = spec1d(x,y,e);

%%
%
% <html><h3>Example 2</h3></html>
%

v.x = x;
v.y = y;
v.e = e;

s_1 = spec1d(v);


%%
%
% <html><h3>Example 3</h3></html>
%

s_3 = spec1d(x,y,e,'x_label','X (points)');
