function [y, p] = mf_interp(xd,yds, x, nb_of_pts, order)
%Lagrange 1D interpolation/smoothing.
%Syntax: [itrp, p] = mf_interp(x,y, {new_x=x, nb_of_pts, order=1}) or mf_interp(y)
%
% This function returns the interpolation/smoothing of (x,y) by Lagrange
% method at 'new_x' values (new_x=x if not precised).
% optional 'nb_of_pts': number of points used for each interpolation computation
% optional 'order': order of polynomial interpolation with nb_of_points
% 'itrp' is the interpolated data. 'p' is last computed polynome.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: Lagrange 1D interpolation/smoothing.

% Part of 'Spectral tools'. E.Farhi. 07/96

% Argument processing -----------------------------------------------
id_tab = [];

if (nargin < 1)
	error('usage: [interpolated,p] = mf_interp(x,y, {new_x=x, nb_of_pts, order=1}) or mf_interp(y)');
end

if (nargin == 1)
	yd=xd;
	xd=1:length(yd);
end

xd = xd(:)'; % make rows
yd = yds(:)';

xdmax = max(abs(xd));
nd=length(xd);

if (nargin <= 4)
	order = [];
end

if isempty(order)
	order = 1;
end

if (nargin <= 3)
	nb_of_pts = [];
end

if isempty(nb_of_pts)
	nb_of_pts = order+1;
end

if (nargin <= 2)
	x = [];
end

if isempty(x)
	x = xd;			% initial X axis
end

x = sort(x);
[xd,X] = sort(xd);
yd = yd(X);
yds = yds(X);

xd = xd  / xdmax;	% xd set to 0:1
x = x(:)';
x = x  / xdmax;

if (length(yd) ~= nd)
	error('x and y vectors must have same number of elements.');
end;
nb_of_pts = min(nb_of_pts, nd);
order = min(order, nb_of_pts-1);
if (nb_of_pts < 2)
	p = [];
	y = yds;
	return;
end
nb_of_pts = nb_of_pts -1;	
fn = floor(nb_of_pts/2);
cn = ceil(nb_of_pts/2);

lx = length(x);

if (lx == nd)
  if(x == xd)			% computing nearest indexes in xd yd to be used
		id_tab = 1:nd;					% no change on x axis.
  else
	dif = x-xd;
	if all( dif == dif(1) )			% x is just xd shifted by constant value
		step = x(2) - x(1);
		id_tab = (1:nd) + round( step/dif(1) );
	end
  end
end
if (isempty(id_tab))				% already computed ?

	for index = 1:lx			% general case.
		[dummy,id] = min(abs(x(index) - xd ));	% closest index in x
		if (length(id))
			id_tab(index) = id(1);
		else
			id_tab(index) = nd;
		end
	end
end

% ----------------------- main iteration loop -----------------------

if (order >= 1)


idprec = 0;
for index = 1:lx			% index in new_x
	id = id_tab(index);		% corresponding index in x data.
	id = max(1, min(id, nd));
	if ((xd(min(id,nd)) == x(index)) & (order == 1) & (nb_of_pts == 2))
		y(index) = yd(id);
else
	if (id ~= idprec)		% need to re-compute polynome.
		idprec = id;
		if (xd(id) < x(index))
			id4 = (-fn):cn;
				% indexes in xd,yd for extracting data
		else
			id4 = (-cn):fn;
		end

		if min(id+id4) <= 1
			id = 1:min(nb_of_pts+1,nd);
		elseif max(id+id4) >= nd
			id = max(1,nd - nb_of_pts):nd;
		else
			id = id+id4;
		end
		xs = xd(id);
		ys = yd(id);

% ------------------- smoothing/interpolation -----------------------

		A = polyfit(xs,ys,order);			% Least Square Approx.

	end;
	
% ----------------------- evaluates polynome -------------------------

	y(index) = polyval(A,x(index));
end
end

else 			% order = 0 -> averaging signal ...
	ys = zeros(1,nd+nb_of_pts);
	xs = ones(1,nd+nb_of_pts)*nb_of_pts;
	for index = 1:nb_of_pts
		ys(index:(nd+index-1)) = ys(index:(nd+index-1)) + yd(1:nd);
		xs(index) = index;
		xs(nd+index) = nb_of_pts-index+1;
	end
	y = ys ./ xs;
	y = y(cn:nd+cn);
	y = y(id_tab);
	A = y(length(y));
end


if (nargout >= 2)	% for standard MatLab/Octave notation.
	p = A;
end


