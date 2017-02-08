function r = minus(s1,s2)
%
% function r = minus(s1,s2)
%
% @SPEC1D/minus function to give the value of spectrum s1 minus a value or spectrum s2.
%
% 1. r = s1-s2     , where s1 and s2 are spectra
% 2. r = s1-a      , where a is to be added to the y axis
% 3. r = s1-[xc,yc], where xc (yc) is added the x (y) axis
%
% If spectra are not the same length, s2 is interpolated to
% the x values of s1.
%
% Minus should not be used to combine counts taken with different
% monitors - use combine(s1,s2)
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

r = plus(s1,s2*(-1));