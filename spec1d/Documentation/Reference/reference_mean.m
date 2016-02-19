%% @spec1d/mean
% This is the reference documentation for the function @spec1d/mean
%
% This function gives the mean in each spectra _s_ where _s_
% can be single or an array. The error _ee_ can also be optionally
% returned. Optional arguments on how to calculate the mean can also be
% given
%%

%% Syntax
%
%    yy = mean(s);
%    [yy, ee] = mean(s);
%    [yy, ee] = mean(s,'method','counts');
%
%% Inputs
%
% * _s_ - Single or vector of spec1d objects
%

%% Optional Inputs
% The optional *method* controls how the mean is calculated. Options are:
% 
% * *mean*: The simple mean of points is calculated.
% * *counts*: The mean is weighted by counts and error.
% * *weight*: The mean is weighted by the error.
% 
% The default is *mean*.
%

%% Outputs
% 
% * _yy_ - Vector of length _s_
% * _ee_ - Vector of length _s_ (optional)

%% Examples
% Theses are examples on using @spec1d/mean
%
% <html><h3>Example 1</h3></html>
%
% Return the mean of each sprecta

s(1) = spec1d(1:5,1:5,0.1*rand(1,5));
s(2) = spec1d(1:5,2:6,0.075*rand(1,5));

yy = mean(s);
fprintf('\ts(1)\ts(2)\n')
fprintf('mean:\t%0.03f\t%0.03f\n',yy)

%%
% 
% <html><h3>Example 2</h3></html>
%
% Return the mean of each sprecta with error
%

[yy, ee] = mean(s);
fprintf('\ts(1)\ts(2)\n')
fprintf('mean:\t%0.03f\t%0.03f\n',yy)
fprintf('error:\t%0.03f\t%0.03f\n',ee)

%%
% 
% <html><h3>Example 3</h3></html>
%
% Return the mean of each sprecta with error by the optional method *weight*
%

[yy, ee] = mean(s,'method','weight');
fprintf('\ts(1)\ts(2)\n')
fprintf('mean:\t%0.03f\t%0.03f\n',yy)
fprintf('error:\t%0.03f\t%0.03f\n',ee)
