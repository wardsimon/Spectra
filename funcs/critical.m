function [y, name, pnames, pin]=critical(x,p, flag)
% critical : Critical
% [y, {name, pnames, pin}]=critical(x,p, {flag}) 
% MFIT Gaussian fitting function
% p = [ dir bkg exp1 ampl1 exp2 ampl2 tc ]

% Author:  CR <c.ruegg@ucl.ac.uk>
% Description: Critical
% set dir=0 for t>tc power law
% fix p(5) and p(6)!
% set dir=1 for t<tc power law
% fix p(5) and p(6)!
% set dir=2 for two power laws

if nargin==2;
    if p(1)==0
        y(abs(x)<p(7))=p(2);
        y(abs(x)>=p(7))=p(2)+p(4)*(x(abs(x)>=p(7))-p(7)).^p(3);
    elseif p(1)==1
        y(abs(x)>p(7))=p(2);
        y(abs(x)<=p(7))=p(2)+p(4)*(p(7)-x(abs(x)<=p(7))).^p(3);
    elseif p(1)==2
        y(abs(x)>p(7))=p(2)+p(6)*(x(abs(x)>p(7))-p(7)).^p(5);
        y(abs(x)<=p(7))=p(2)+p(4)*(p(7)-x(abs(x)<=p(7))).^p(3);
    end
    y=y';
else
	y=[];
	name='critical';
	pnames=str2mat('dir','bkg','exp1','ampl1','exp2','ampl2','tc');
	if flag==1, pin=[0 0 0 0 0 0 0]; else pin = p; end
	if flag==2
	end
end
