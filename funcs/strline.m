function [y, name, pnames, pin]=strline(x,p, flag)
% strline   : slope/line
% function [y, {name, pnames, pin}]=strline(x,p, {flag})
%
% MFIT slope/line function ax+b
% p = [ a b ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description:  slope/line

if nargin==2;
	y=p(2)+p(1)*x;
else
	y=[];
	name='Straight line';
	pnames=str2mat('Gradient','Background');
	if flag==1, pin=[0 0]; else pin = p; end
	if flag==2
    	mf_msg('Click on point 1');
      [x1 y1]=ginput(1);
		mf_msg('Click on point 2');
      [x2 y2]=ginput(1);
      grad=(y2-y1)/(x2-x1);
      bg=y1-(x1*grad);
      pin=[grad bg];
    end
end
