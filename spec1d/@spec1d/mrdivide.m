function r=mrdivide(s1,s2)
%
% function r=mrdivide(s1,s2)
%
% @SPEC1D/RDIVIDE Rdivide for 1D spectra
%
% Valid expressions: 
%
% 1. r=s1./s2      , where s1 and s2 are spectra
% 2. r=s1./yc      , where yc is a constant that divides the y-axis
% 3. r=s1./[xc,yc] , where xc (yc) divides the x (y) axis by a constant
%
% DFM 1.4.98
%
r=rdivide(s1,s2);