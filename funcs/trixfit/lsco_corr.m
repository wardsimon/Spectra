function [yout]=trixfit_corr(x,y)
%
% TRIXFIT function to apply a correction factor to the
%         calculated intensity
%
% Des McMorrow, June 2001
%

%----- Correct for lambda/2 for ki scans. 
%      This needs to be checked for different instruments
%
%yout=y.*(1+91.34./((x+5).^3.58));

yout=y;
