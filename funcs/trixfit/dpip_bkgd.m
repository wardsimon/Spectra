function [bcgd]=trixfit_bkgd(x,p)
%
% TRIXFIT function to define the background
% 
% Des McMorrow, June 2001
% Beni Thielemann, October 2006
bcgd=p(15)+p(16)*x.^2+p(17)*exp(-p(18)*x);
