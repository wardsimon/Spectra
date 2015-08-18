function [y, name, pnames, pin] = green(x,p, flag)
% green     : Green function
% [y, name, pnames, pin] = green(x,p, flag)
%
% Green function (two symetrical peaks)
%    p = [ intensity doublet_pos full_width bg ]
% maximum is I/full_width at x=doublet_pos
% and at x=0, y=I*full_width/(x=doublet_pos)^2
% integral is pi*p(1) when bg = 0

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  Green function

if nargin==2	% Bose becomes 1/omega for a given T
	y = (abs(p(3)*p(1))*p(2)^2 ) ./ ( (p(2)^2 - x.^2).^2 + (x*p(3)).^2) + p(4);
else
	y=[];
	name='Green';
	pnames=str2mat('Intensity','Doublet.pos', 'Width','Background');
	if flag==1, pin=[1 0.1 1 0]; else pin = p; end
	if flag==2
		mf_msg('Click on background');
		[cen bg]=ginput(1);
		mf_msg('Click on doublet');
		[cen amp]=ginput(1);
		mf_msg('Click on doublet width');
		[width y]=ginput(1);
		width=abs(width-cen);
		pin=[amp cen width bg ];
	end

end



