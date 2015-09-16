function [y, name, pnames, pin]=quadrat(x,p, flag)
% quadrat   : quadratic
% function [y, {name, pnames, pin}]=quadrat(x,p, {flag})
%
% MFIT quadratic function a*x^2+b*x+c
% p = [ a b c ]

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  quadratic

if nargin==2;
	y=p(3)+p(2)*x+p(1)*x.^2;
else
	y=[];
	name='Quadratic ax^2+bx+c';
	pnames=str2mat('a','b','c');
	if flag==1, pin=[1 1 1]; else pin = p; end
	if flag==2
		mf_msg('Click on point 1');
		[x1 y1]=ginput(1);
		mf_msg('Click on point 2');
		[x2 y2]=ginput(1);
		mf_msg('Click on point 3');
		[x3 y3]=ginput(1);
		[pin,ny] = polyfit([x1 x2 x3],[y1 y2 y3],2);
%		a=po(1); b=po(2); c=po(3);
    end
end
