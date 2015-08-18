function [y, name, pnames, pin]=voigt(x,p, flag)
% voigt     : Voigt
% function [y, {name, pnames, pin}]=voigt(x,p, {flag})
%
% MFIT Voigt fitting function
% p = [ Amplitude Centre Gauss_width Lorz_width Background ]

% Author:  MZ <mzinkin@sghms.ac.uk> adapted from DFM
% Description:  Voigt

if nargin==2;
	N = 16;
	b = -sqrt(log(2))/p(3);
	a = b*p(4);
	b = b*2*1i;
	z = a + b*(x-p(2));

	M=2*N; M2=2*M; k=[-M+1:1:M-1]';
	L=sqrt(N/sqrt(2));
	tt=(L*tan(k*pi/M2)).^2;
	f=[0; exp(-tt).*(L^2+tt)];
	a=real(fft(fftshift(f)))/M2;
	a=flipud(a(2:N+1));
	l=L-z;
	Z=(L+z)./l;
	pp=polyval(a,Z);
	y=p(5)+p(1)*real(2*pp ./l.^2+(1/sqrt(pi))*ones(size(z)) ./l);
else
	y=[];
	name='Voigt function';
	pnames=str2mat('Amplitude ','Centre','Gauss._width',...
						'Lorz._width','Background');
	if flag==1, pin=[0 0 1 1 0]; else pin = p; end
	if flag==2
		mf_msg('Click on peak');
		[cen amp]=ginput(1);
		mf_msg('Click on width');
		[width y]=ginput(1);
		width=abs(width-cen);
		mf_msg('Click on background');
		[x bg]=ginput(1);
		amp=amp-bg;
		pin=[2*amp cen width/2 width/2  bg];
	end
end
