function [y, name, pnames, pin]=gaussx2(x,p, flag)
% gaussx2   : 2 gaussians
% [y, {name, pnames, pin}]=gaussx2(x,p, {flag})
%
% MFIT 2 gaussians fitting function. Centres are relative to
%	centre of first gaussian
% p = [ Amp1 Centre1 Width1 Amp2 Centre2 Width2 BakG ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description:  2 gaussians

global ngss;

if nargin==2
    y=p(1)*exp(-0.5*(((x-p(2))/p(3)).^2)) ...
      + p(4)*exp(-0.5*(((x-p(5))/p(6)).^2)) ...
      +p(7);
else
	y=[];
	if flag==1, pin=[0 0 1 0 0 1 0]; else pin = p; end
	name='2 Gaussians';
	pnames=str2mat('Amplitude_1','Centre_1','Width_1',...
						'Amplitude_2','Centre_2','Width_2',...
						'Background');

   if flag==2

      mf_msg('Click on background');    
		[x bg]=ginput(1);

		mf_msg('Click on peak 1');
		[cen amp]=ginput(1);
		mf_msg('Click on width 1');
		[width y but]=ginput(1);
		width=abs(width-cen);
		amp=amp-bg;
		pin=[amp cen width];

		mf_msg('Click on peak 2');
		[cen amp]=ginput(1);
		mf_msg('Click on width 2');
		[width y but]=ginput(1);
		width=abs(width-cen);
		amp=amp-bg;
		pin=[pin amp cen width];

		pin=[pin bg];
	end

end
