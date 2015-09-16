function [y, name, pnames, pin] = sdk(x,p,flag)
% sdk       : 2 Airy functions product power 3 for SDK
% [y, name, pnames, pin] = sdk(x,p,flag)
% % 2 Airy functions product for Fabry-Perot with
% p=[ intensity M1 Asym1 ISL1 M2 Asym2 ISL2 Xelast BckG ]

% EF 09.97

if (nargin == 2)
	delt1=-p(3)-4*p(2)*p(3)^3;
	delt2=-p(6)-4*p(5)*p(6)^3;

	phi1=pi*((x-p(8)+delt1)/p(4));
	phi2=pi*((x-p(8)+delt2)/p(7));


	d1=(1+p(2)*(sin(phi1)).^2);
	i1=((1./d1)-p(3)*(p(2)*sin(2*phi1))./d1.^2);

	d2=(1+p(5)*(sin(phi2)).^2);
	i2=((1./d2)-p(6)*(p(5)*sin(2*phi2))./d2.^2);

	y=(p(1).*i1.*i2).^3+p(9);
else
	y=[];
	name='2 Airy functions product for SDK';
	pnames=str2mat('Intensity','M1(width)', 'Asym1', 'Period1(ISL)','M2(width)', 'Asym2', 'Period2(ISL)','Centre','Background');
	if flag==1, pin=[10 150 0.01 270 150 0.01 270 512 1 ]; else pin = p; end
	if flag==2
		mf_msg('Click on central peak');
		[cen amp]=ginput(1);
		mf_msg('Click on typical width');
		[width1 y]=ginput(1);
		width1=abs(width1-cen);
		mf_msg('Click on Ghost1.');
		[ISL1 amp]=ginput(1);
		ISL1 = abs(ISL1 - cen);
		mf_msg('Click on Ghost2.');
		[ISL2 y]=ginput(1);
		ISL2 = abs(ISL2 - cen);
		mf_msg('Click on Background.');
		[dummy bg]=ginput(1);
		pin=[ amp^(1/3) width1 0 ISL1 width1 0 ISL2 cen bg];
	end
end
