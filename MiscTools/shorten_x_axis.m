function shr_x = shorten_x_axis (x, len)
% shorten_x_axis : Shortens vector's ends.
%Syntax: shr_x = extend_x_axis (x, len)
%
% Shortens the x axis by 'len' points, distributed symetricaly on
% both sides, assumed step is constant.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: remove some data at vector's ends.

% uses : none
% E.Farhi 12/95

if (nargin < 2)
	error('usage : shr_x = shorten_x_axis (x, len)');
end

lx = length(x);
fn = floor(len/2);

shr_x = x((fn+1):(lx-len+fn));

%end



