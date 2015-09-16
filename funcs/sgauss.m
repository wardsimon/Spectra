function [y, name, pnames, pin]=sgauss(x,p, flag)
% sgauss    : Gaussian plus a sloping background
% function [y, name, pnames, pin]=sgauss(x,p, flag)
%
% MFIT Gaussian plus a sloping background function
% p = [ Constant_b/g Slope_b/g Amplitude Centre Width ]

% Author:  AS 12.12.94
% Description:  Gaussian + sloping

if nargin==2;
	y=p(1)+x*p(2)+p(3)*exp(-0.5*((x-p(4))/p(5)).^2);
else
	y=[];
	name='Gaussian + sloping background';
	pnames=str2mat('Background','Slope.bkg','Amplitude','Centre','Width');
	if flag==1, pin=[1 0 0 0 1]; else pin = p; end
	if flag==2
		mf_msg('Click on first point on bg');
		[xi1 yi1]=ginput(1);
		mf_msg('Click on second point on bg');
		[xi2 yi2]=ginput(1);
		m=(yi2-yi1)/(xi2-xi1);
		c=yi1-m*xi1;
		mf_msg('Click on peak');
		[cen amp]=ginput(1);
		amp=amp-(m*cen+c);
		mf_msg('Click on width');
		[width y]=ginput(1);
		width=abs(width-cen);
		pin=[c m amp cen width];
	end
end
