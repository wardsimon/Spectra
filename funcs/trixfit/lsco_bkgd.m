function [bcgd]=trixfit_bkgd(x,p)
%
% TRIXFIT function to define the background
% 
% Des McMorrow, June 2001
%
bcgd=p(11)+p(12)*x.^2+p(13)*exp(-p(14)*x);
