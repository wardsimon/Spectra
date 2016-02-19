%% @spec1d/extract
%
% This function extracts raw data from spec1d object _s_, or from
% spec1d array. If all spectra in _s_  have same number of data points a
% matrix is returned, else a cell of length _s_.
%%

%% Syntax
%
%    [x, y, e] = extract(s)

%% Inputs
%
% * _s_ - Single or vector of spec1d objects
%

%% Outputs
% 
% <html><ul style="list-style-type:square">
%  <li><i>x</i> - X-points</li>
%  <ul> 
%  <li><i>vector</i> - Column vector.</li>
%  <li><i>matrix</i> - MxN matrix where M is number of points and N is the length of _s_.</li>
%  <li><i>cell</i> - Cell of column vecors, with the length of _s_.</li></ul>
%  <li><i>y</i> - Y-points</li>
%  <ul> 
%  <li><i>vector</i> - Column vector.</li>
%  <li><i>matrix</i> - MxN matrix where M is number of points and N is the length of _s_.</li>
%  <li><i>cell</i> - Cell of column vecors, with the length of _s_.</li></ul>
%  <li><i>e</i> - Error associated to Y-points</li>
%  <ul> 
%  <li><i>vector</i> - Column vector.</li>
%  <li><i>matrix</i> - MxN matrix where M is number of points and N is the length of _s_.</li>
%  <li><i>cell</i> - Cell of column vecors, with the length of _s_.</li></ul>
% </ul>
% </html>
% 

%% Example
% This is an example on using @spec1d/extract
%
% <html><h3>Example 1</h3></html>
%

s = spec1d(1:0.2:5,sin(linspace(0,2*pi,21)),rand(21,1));
[x, y, e] = extract(s);

fprintf('x\ty\te\n')
fprintf('%0.03f\t%0.03f\t%0.03f\n',[x,y,e]')

%% See Also
% <html><a href="{{ site.url }}/@spec1d.get/index.html">Get</a></html>