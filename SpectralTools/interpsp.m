function [y, p] = interpsp(xd,yds, x, nb_of_pts, order, last)
% interpsp : Lagrange 1D interpolation/smoothing.
%Syntax: [itrp, p] = interpsp(x,y, {new_x=x, nb_of_pts, order=1,mode='silent'}) or interpsp(y)
%
% This function returns the interpolation/smoothing of (x,y) by Lagrange
% method at 'new_x' values (new_x=x if not precised).
% optional 'nb_of_pts': number of points used for each interpolation computation
% optional 'order': order of polynomial interpolation with nb_of_points
% 'itrp' is the interpolated data. 'p' is last computed polynome.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: Lagrange 1D interpolation/smoothing.

% Part of 'Spectral tools'. E.Farhi. 07/96 rev 11/97
% uses : 

% Argument processing -----------------------------------------------

if ~exist('last')
	last = 'silent';
end

if isempty(last)
	last = 'silent';
end

if strcmp(last,'silent')
	tmp = 0;
else
	tmp = 1;
end

if (tmp>0)
	fprintf(1,'Smoothing/Interpolation\n');
end

if (nargin < 1)
	error('usage: [interpolated,p] = interpsp(x,y, {new_x=x, nb_of_pts, order=1}) or interp(y)');
end

if (nargin == 1)
	yd=xd;
	xd=1:length(yd);
end

xd = xd(:);
yd = yds(:);

xdmax = max(abs(xd));
nd=length(xd);

if ~exist('order')
	order = [];
end
if isempty(order)
	order = 1;
end

if ~exist('nb_of_pts')
	nb_of_pts = [];
end
if isempty(nb_of_pts)
	nb_of_pts = order+1;
end

if ~exist('x') 
	x = [];
end
if isempty(x)
	x = xd;			% initial X axis
	if (tmp>0)
		disp('new_x = x');
	end
end

if (xd(1) > xd(nd))
	xd = xd(nd:-1:1);
	yd = yd(nd:-1:1);
	x = x(length(x):-1:1);
	xddec = 1;
else
	xddec = 0;
end

xd = xd  / xdmax;	% xd set to 0:1
x = x(:);
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
	
if (tmp>0)
	fprintf(1,'Nb of points %i, Order %i. ',nb_of_pts,order);
end
nb_of_pts = nb_of_pts -1;
fn = floor(nb_of_pts/2);
cn = ceil(nb_of_pts/2);

lx = length(x);

if (tmp>0)
	fprintf(1,'Running ');
end
t0 = clock;
X = eye(order+1);
Y = zeros(1,order+1);

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
if (~exist('id_tab'))				% already computed ?
	x = sort(x);
	xd = sort(xd);
	yds = sort(yds);
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

if (tmp>0)
	fprintf(1,'(10 %%) #');
end

idprec = 0;
for index = 1:lx			% index in new_x
	id = id_tab(index);		% corresponding index in x data.
	id = max(1, min(id, nd));

	if (id ~= idprec)		% need to re-compute polynome.
		idprec = id;
		if (xd(id) < x(index))
			id = max(1, id-fn):min(nd,id+cn);	% indexes in xd,yd for extracting data
		else
			id = max(1, id-cn):min(nd,id+fn);
		end
		if id == 1
			id = 1:min(nb_of_pts+1,nd);
		elseif id == nd
			id = max(1,nd - nb_of_pts):nd;
		end
		xs = xd(id);
		ys = yd(id);

% ------------------- smoothing/interpolation -----------------------

		for j=1:(order+1)
			X(j,j) = sum(xs.^(2*j-2));
			Y(j) = sum(ys.*(xs.^(j-1)));
			for k=(j+1):(order+1)
				X(k,j) = sum(xs.^(k+j-2));
				X(j,k) = X(k,j);
			end;
		end;
		
		A = X\Y';			% Least Square Approx.

	end;
	
% ----------------------- evaluates polynome -------------------------

	k=0;
	for j=1:order+1
		k = k +A(j)*x(index)^(j-1);
	end;
	y(index) = k;

	if (rem(index,ceil(lx/10)) == 0)
		if (tmp>0)
			fprintf(1,'#');
		end
	end
end


else 			% order = 0 -> averaging signal ...
	ys = zeros(1,nd+nb_of_pts);
	xs = ones(1,nd+nb_of_pts)*nb_of_pts;
	for index = 1:nb_of_pts
		ys(index:(nd+index-1)) = ys(index:(nd+index-1)) + yd(1:nd);
		xs(index) = index;
		xs(nd+index) = nb_of_pts-index+1;
		if (tmp>0)
			fprintf(1,'#');
		end
	end
	y = ys ./ xs;
	y = y(cn:nd+cn);
	y = y(id_tab);
	A = y(length(y));
end

if (tmp>0)
	fprintf(1,' - OK %5.1f s.\n',etime (clock, t0));
end

if (nargout >= 2)	% for standard MatLab/Octave notation.
	X = 1 / xdmax;
	for index = 1:order+1
		A(index) = A(index) * X^(index-1);
	end
        p = A(:);
	p = fliplr(p);
	if (tmp>0)
		disp('Interpolation polynome (Last) :');
		disp(p);
	end;
end

if (size(yds,2)==1)
	y = y';
end

if (xddec == 1)
	y = y(length(y):-1:1);
end
