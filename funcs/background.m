function [y, name, pnames, pin]=background(x,p, flag)
% background: constant
% function [y, {name, pnames, pin}]= background(x,p, {flag})
%
% MFIT constant

% Author:  EF
% Description:  slope/line

if nargin==2;
	y=p(1)*ones(size(x));
else
	y=[];
	name='Constant Background';
	pnames=str2mat('Background');
	if flag==1, pin=1; else pin = p; end
	if flag==2
		mf_msg('Click on background');
		[x1 y1]=ginput(1);
		pin=y1;
	end
end
