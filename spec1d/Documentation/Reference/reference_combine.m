%% @spec1d/combine
% This is the reference documentation for the function @spec1d/combine.
%
% This function to combines two or more spectra. If the x values of two points differ by less
% than tolerance _toll_, then the points are combined.
%%

%% Syntax 
% 
%   s_out = combine(toll,s,varargin)
%

%% Inputs
% 
% * _toll_ - single value
% * _s_ - single spec1d object or array
% * _varargin_ - Optional Parameter/Value pairs
% 

%% Optional Inputs
% 
% Depending on the optional *method*, points are combined as
%
% # *mean*    : Simple means for x and y, errors are averaged in quadrature.
% # *counts*	: Restablishes normalisation and original counts assuming
%             square-root statistics, is correct for normalised counts
% # *weight*	: Weights to inverse error. For more general data.
%
% Default is *counts*
%
% Depending on the optional *indexing*, points are indexed as
%
% # *relative*    : Histogram bining of size 'toll'. This is the same as
%                 using the rebin function except for multiple spectum.
% # *absolute*	: The gap between the points is  always greater then 'toll'
%                 before any averaging.
% # *legacy*      : Replicate the original @spec1d/combine
%
% Default is *relative* due to speed considerations.

%% Outputs
% 
% * _s_out_ - spec1d object

%% Examples
% These are some examples on using @spec1d/combine. Combining spectra _s1_, _s2_ and _s3_ when x values differ by less than 0.5.
% 
% <html><h3>Example 1</h3></html>
%
% Combine _s1_, _s2_, _s3_ by counts method
%

s1 = spec1d(1:0.2:5,sin(linspace(0,2*pi,21)),0.1);
s2 = spec1d((1:0.2:5) + 0.2*rand(1,21),(sin(linspace(0,2*pi,21))) - 0.2*rand(1,21),0.1);
s3 = spec1d((1:0.2:5) - 0.2*rand(1,21),(sin(linspace(0,2*pi,21))) + 0.2*rand(1,21),0.1);

s_out1 = combine(0.5,s1,s2,s3);

figure
plot(s1,s2,s3,s_out1)

%% 
%
% <html><h3>Example 2</h3></html>
%
% Combine by the 'mean' method
%

s_out2 = combine(0.5,s1,s2,s3,'method','mean');

figure
plot(s1,s2,s3,s_out2)

%% 
% 
% <html><h3>Example 3</h3></html>
%
% Combine _s1_, _s2_, _s3_ by weight method and indexing 'absolute'
%

s_out3 = combine(0.5,s1,s2,s3,'method','weight','indexing','absolute');

figure
plot(s1,s2,s3,s_out3)

%% Notes
% The binning method and indexing results in different final objects. It is
% important to note the advantages/disadvantages of each method. 

figure
hold on
plot(get(s1,'x'),get(s1,'y'),'r-')
plot(s_out1,s_out2,s_out3)

%% See Also
% <html><a href="{{ site.url }}/@spec1d.rebin/index.html">Rebin</a></html>