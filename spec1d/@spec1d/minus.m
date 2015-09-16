function r=minus(s1,s2)
%
% function r=minus(s1,s2)
%
% @SPEC1D/MINUS Minus for 1D spectra
%
% Valid expressions: 
%
% 1. r=s1-s2     , where s1 and s2 are spectra
% 2. r=s1-a      , where a is to be subtracted from the y axis 
% 3. r=s1-[xc,yc], where xc (yc) is subtracted from the x (y) axis 
%
% DFM 1.4.98
r=plus(s1,s2*(-1));