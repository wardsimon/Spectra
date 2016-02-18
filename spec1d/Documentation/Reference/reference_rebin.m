%% @spec1d/rebin
% This is the reference documentation for the function @spec1d/rebin.
%
% This function re-bins a spec1d spectrum _s_ to the bins specified by _dx_.
%%

%% Syntax 
% 
%   s_out = rebin(s,dx,varargin)
%

%% Inputs
% 
% <html><ul style="list-style-type:square">
%  <li><i>s</i> - Single spec1d object or array</li>
%  <li><i>dx</i> - Binning argument</li>
%  <ul> 
%  <li><i>scalar</i> - Specifies bin widths.</li>
%  <li><i>vector</i> - Specifies the bin centers.</li>
%  <li><i>sped1d</i> - Uses x-values of dx as bin centers</li>
%  <li><i>string</i> - Use a binning algorithm described in the documentation histcounts.</li></ul>
%  <li><i>varargin</i> - Optional Parameter/Value pairs</li>
% </ul>
% </html>
% 

%% Optional Inputs
% The optional *method* controls how the x-value of the bin is interpreted.
% Options are:
%
% # *average* : New x values are average within bins (weighted by error)
% # *force*   : New x values are as specified in dx
% # *interp*  : New x-values as specified by dx, y and e interpolated from 'average'-x,
%
% The default is *average*
%
% The methods above can also have the *byError* flag.
%                    
% # *0* : Data is binned by x. 'Default'
% # *1* : Data is binned to a constant error. Only *method*, *average* is supported.

%% Outputs
% 
% * _s_out_ - spec1d object

%% Examples
% These are some examples on using @spec1d/rebin. Rebin spectra _s1_.
%
% <html><h3>Example 1</h3></html>
%
% Rebin _s_ to x-value steps 0.5 
%

s = spec1d(1:0.2:5,sin(linspace(0,2*pi,21)),0.1);
s_out1 = rebin(s,0.5);

figure
plot(s,s_out1)

%% 
%
% <html><h3>Example 2</h3></html>
%
% Rebin _s_ to an x vector
%

x = 1:0.4:5;
x(3:4:end) = []; 
s_out2 = rebin(s,x);

figure
plot(s,s_out2)

%% Notes
% The width of boundary bins is twice the distance to neighbour bin.

%% See Also
% <html><a href="{{ site.url }}/@spec1d.combine/index.html">Combine</a></html>