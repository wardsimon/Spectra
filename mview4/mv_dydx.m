function [dydx, xd, err]=mv_dydx(x, y, yerr)
%----------------------------------------------------------------
%function [dydx, xd, err]=deriv(x, y, yerr)
%  Returns derivative of vector y wrt x
%
% M. Zinkin 12.12.94
%----------------------------------------------------------------

ym=[y(:); NaN; NaN];
yp=[NaN; NaN; y(:)];
y =[NaN; y(:); NaN];

xm=[x(:); NaN; NaN];
xp=[NaN; NaN; x(:)];
x =[NaN; x(:); NaN];

em=[yerr(:); NaN; NaN];
ep=[NaN; NaN; yerr(:)];
e =[NaN; yerr(:); NaN];

dydx= 0.5*((yp-y)./(xp-x) + (y-ym)./(x-xm));
err = 0.25*( (ep.^2+e.^2)./(xp-x).^2 + ...
             (e.^2+em.^2)./(x-xm).^2);
xd=(xp+xm+2*x)/4;

l=length(dydx)-2;
dydx=dydx(3:l);
err=err(3:l);
xd=xd(3:l);
