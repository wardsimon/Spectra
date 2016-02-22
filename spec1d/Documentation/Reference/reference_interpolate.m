%% @spec1d/interpolate
% This is the reference documentation for the function @spec1d/interpolate
%
% This function is to interpolate the given spectra _s_ for the x-ranges given by the vector _x_new_.
%%

%% Syntax
%
%    s_out = interpolate(s,x_new,options)

%% Inputs
%
% <html><ul style="list-style-type:square">
%  <li><i>s</i> - Single spec1d object or array</li>
%  <li><i>x_new</i> - New x-axis</li>
%  <ul> 
%  <li><i>vector</i> - Specifies the new x-points.</li>
%  <li><i>sped1d</i> - Use x-values from the spec1d object as the new x-points.</li>
% </ul>
% </html>
%

%% Optional Inputs
%
% Depending on the optional *method*, points are interpolated as
% * *weightedpoly* : Adaptive polynomial interpolation where errors are calculated and minimised for a n-length polynomial function.
% * *linear*	   : Nearest neighbout interpolation where values are weighted by inverse distance (including errors)
% * *builtin*	   : Use the MATLAB interp1 function. Optional arguments can be passed
% * *FUNCTION*     : The value FUNCTION is a function which can be explained with the documentation "help('sdinterp')". Optional arguments can be passed.
%
% Default is *weightedpoly*. Use 'builtin' or 'linear' for speed or non-uniform data.
%
% If the *method* *weightedpoly* is selected, the optional *order*, parameter
% is available. The value is a positive integer and corresponds to the
% order of the fitted polynomial. Default is -1, which optimises for least
% error.
%

%% Outputs
%
% * _s_out_ - Spec1d object of the same size as the input array

%% Notes
%
% * _x_new_ is a vecor of arbitary length where min(_x_new_) < min(_s_.x) and
% max(_x_new_) < max(_s_.x)
% * Error calculation is accurate when square-root statistics is valid. Otherwise the error is an approximation. 

%% Examples
% These are  examples using @spec1d/interpolate
%
% <html><h3>Example 1</h3></html>
%

s = spec1d(linspace(0,2*pi,10),sin(linspace(0,2*pi,10)),0.1);
x_new = 0.1:0.1:(2*pi);
s_1 = interpolate(s,x_new);

figure
plot(s);
hold on
plot(s_1,'LineStyle','-','MarkerSize',4)
legend({'s','interpolate(s)'})

%%
%
% <html><h3>Example 2</h3></html>
%
% Specify the order of the polynominal to be used.
%

s_2 = interpolate(s,x_new,'order',3);

figure
plot(s);
hold on
plot(s_2,'LineStyle','-','MarkerSize',4)
legend({'s','interpolate(s)'})

%%
%
% <html><h3>Example 3</h3></html>
%
% Use the linear interpolation method.
%

s_3 = interpolate(s,x_new,'method','linear');

figure
plot(s);
hold on
plot(s_3,'LineStyle','-','MarkerSize',4)
legend({'s','interpolate(s)'})

%%
%
% <html><h3>Example 5</h3></html>
%
% Interpolation using a function in the sdinterp library
%


s_4 = interpolate(s,x_new,'method','cbezier');

figure
plot(s);
hold on
plot(s_4,'LineStyle','-','MarkerSize',4)
legend({'s','interpolate(s)'})