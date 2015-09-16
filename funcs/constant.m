function [y, name, pnames, pin]=constant(x,p, flag)
% constant  : constant
% function [y, {name, pnames, pin}]= constant(x,p, {flag})
%
% MFIT constant

% Author:  EF
% Description:  slope/line

if nargin==2;
	y=p(1)*ones(size(x));
else
	y=[];
	name='Constant';
	pnames=str2mat('Constant');
	if flag==1, pin=1; else pin = p; end
	if flag==2
		mf_msg('Click on constant(y axis)');
		[x1 y1]=ginput(1);
		pin=y1;
	end
end
