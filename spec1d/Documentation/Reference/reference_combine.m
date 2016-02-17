%% @spec1d/combine
% This is a the reference documentation for the function @spec1d/combine
% This function to combines two or more spectra. If the x values of two points differ by less
% than tolerance 'toll', then the points are combined.
%%

%% Usage 
% 
%  # s_out = combine(toll,s,varargin)
%

%% Inputs
% 
% * _toll_ - single value
% * _s_ - single spec1d object or array
% * _varargin_ - Optional Parameter/Value pairs
% 

%% Optional Inputs
% 
% Depending on the optional 'method', points are combined as
% 'mean'    : Simple means for x and y, errors are averaged in quadrature.
% 'counts'	: Restablishes normalisation and original counts assuming
%             square-root statistics, is correct for normalised counts
% 'weight'	: Weights to inverse error. For more general data.
% Default is 'counts'
%
% Depending on the optional 'indexing', points are indexed as
% 'relative'    : Histogram bining of size 'toll'. This is the same as
%                 using the rebin function except for multiple spectum.
% 'absolute'	: The gap between the points is  always greater then 'toll'
%                 before any averaging.
% 'legacy'      : Replicate the original @spec1d/combine

%% Outputs
% 
% * _s_out_ - spec1d object

%% Notes
% 
% * Default 'indexing' is 'relative' due to speed considerations.
% 

%% Examples
% These are some examples on using @spec1d/combine. Combining spectra s1,s2 and s3 when x values differ by less than 0.5.

%%% Example 1
% Combine s1,s2,s3 by counts method
s1 = spec1d(1:0.2:5,sin(linspace(0,2*pi,21)),0.1);
s2 = spec1d((1:0.2:5) + 0.2*rand(1,21),(sin(linspace(0,2*pi,21))) - 0.2*rand(1,21),0.1);
s3 = spec1d((1:0.2:5) - 0.2*rand(1,21),(sin(linspace(0,2*pi,21))) + 0.2*rand(1,21),0.1);

s_out = combine(0.5,s1,s2,s3);

figure
plot(s1,s2,s3,s_out)

% 
% > % 
% >r = combine(0.01,s1,s2,s3,'method','mean') % combine s1,s2,s3 by mean method
% >r = combine(0.01,s1,s2,s3,'method','weight','indexing',absolute) % combine s1,s2,s3 by
% weight method and indexing 'absolute'