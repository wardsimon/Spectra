function [bragg]=rc_bragg(M)
%
% MATLAB function to calculate the widths (FWHM)
% of a Bragg peak from the resolution matrix M
%
% Output: bragg, Qx, Qy, Qz, Vanadium and DEE widths
% 
% A.T.
%
bragg(1)=2.3548/sqrt(M(1,1));
bragg(2)=2.3548/sqrt(M(2,2));
bragg(3)=2.3548/sqrt(M(3,3));
[r,bragg(4)]=rc_phon(1,M,[0 0 0 1]);
bragg(5)=2.3548/sqrt(M(4,4));

