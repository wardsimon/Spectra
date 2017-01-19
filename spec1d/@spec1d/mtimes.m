function r=mtimes(s1,s2)
%
% function r=times(s1,s2)
%
% @SPEC1D/TIMES Times for 1D spectra
%
% Valid expressions: 
%
% 1. r=s1.*s2      , where s1 and s2 are spectra
% 2. r=yc.*s1      , where yc multiplies the y axis
% 3. r=[xc,yc].*s1 , where xc (yc) multiplies the x (y) axis
%
r=times(s1,s2);