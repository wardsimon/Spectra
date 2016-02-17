%% @spec1d/cut 
% This is a the reference documentation for the function @spec1d/cut
%
% This function Cuts data from a spec1d spectrum or spec1d array _s1_ using the 
% xranges specified by [xleft1 xright1], etc.
%%

%% Syntax 
% Cut data outside of the range xleft-xright
%
%    s_out = cut(s,[xleft xright]); 
%    s_out = cut(s,[xleft1 xright1],[xleft2 xright2],....); 
%
% Cut data in the range xleft-xright
%
%    s_out = cut(s,[xright xleft]);

%% Inputs
% 
% * _s_ - single or array of spec1d objects
% * _[xleft xright]_ cutting ranges. Multpile ranges can be specified one
% afer another.
% 

%% Outputs
% 
% * _s_out_ - spec1d object of thge same size as s

%% Notes
% 
% * If spectra are not the same length, _s2_ is interpolated to the x
% values of _s1_
% 

%% Examples
% These are some examples on using @spec1d/cut

%% Example 1
% Cut everything before an x-value of 5; 

s = spec1d(1:10,rand(10,1),0.1);
s_1 = cut(s,[NaN 5]);
figure
plot(s,s_1)
legend({'s','s_1'})

%% Example 2
% Cut everything after an x-value of 5; 
s_2 = cut(s,[5 NaN]);
figure
plot(s,s_2)
legend({'s','s_2'})
snapnow

%% Example 3
% Cut out the range 2.5 to 7.5
s_3 = cut(s,[7.5 2.5]);
figure
plot(s,s_3)
legend({'s','s_3'})

%% Example 4
% Leave only in the range 2.5 to 7.5
s_4 = cut(s,[2.5 7.5]);
figure
plot(s,s_4)
legend({'s','s_4'})