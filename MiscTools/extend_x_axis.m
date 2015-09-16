function ext_x = extend_x_axis (x, len)
% extend_x_axis : Extend vector at ends
%Syntax: ext_x = extend_x_axis (x, len)
% Extends the x axis by 'len' points, distributed symetricaly on
% both sides, assumed step is constant.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: add some data at vector's ends.

% uses : vect2row.m
% E.Farhi 12/95

if (nargin < 2)
	error('usage : ext_x = extend_x_axis (x, len)');
end

xr = vect2row(x);

lx = length(x);
fn = floor(len/2);
if (lx == 1)
	step =1;
else
	step = (x(lx) - x(1)) / (lx - 1);
end

i1 = x(1) - step * fn;
i4 = x(lx) + step * (len-fn);

ext_x = linspace(i1, i4, lx+len);

if (is_column(x))
	ext_x = ext_x';
end
%end



