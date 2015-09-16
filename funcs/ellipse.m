function [y, name, pnames, pin]=ellipse(x,p, flag)
% ellipse   : Ellipse
% [y, {name, pnames, pin}]=ellipse(x,p, {flag})
%
% MFIT Ellipse  fitting function
%	MZ 29.11.94
% p = [ a b theta(deg) x0 y0 ]
%
if nargin==2
	a=1/p(1)^2;
	c=1/p(2)^2;
	t=p(3)*pi/180;
	x0=p(4);
	y0=p(5);

	u=a*cos(t)^2+c*sin(t)^2;
	v=0.5*(a-c)*sin(2*t);
	w=a*sin(t)^2+c*cos(t)^2;
	x=x-x0;

	y=(-v*x+sqrt(v^2*x.^2-w*u*x.^2+w))/w+y0;
else
	y=[];
	name='ellipse';
	pnames=str2mat('Pr. axis a','Pr. axis b','Theta','x0','y0');
	if flag==1, pin=[1 1 0 0 0]; else pin = p; end
	if flag==2
		mf_msg('Click on centre');
		[x0 y0]=ginput(1);
		mf_msg('Click on Pr. axis 1');
		[x1 y1]=ginput(1);
		mf_msg('Click on Pr. axis 2');
		[x2 y2]=ginput(1);
		pin=[sqrt((x1-x0)^2+(y1-y0)^2) sqrt((x2-x0)^2+(y2-y0)^2) ...
                     180/pi*atan((y1-y0)/(x1-x0)) x0 y0];
	end
end
