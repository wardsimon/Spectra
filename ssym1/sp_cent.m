function [centre, cencriteria,cencrit,dc] = sp_cent(x,y,peak)
% [centre, cencriteria,cencrit] = sp_cent(x,y,peak);
% search symetry in spectrum

% gravity centre : sum(y.*x)./sum(y)
% other method : for xo, sum(y(x)-y(2*xo-x))^2 to be minimized.

% first we need to evaluate the fluctuations (noise) in y.
x=x(:); y=y(:);
ly = length(y);
ny = median(abs(y-smooth(y)));
dx = max(x)-min(x);
index = find( (peak(:,5) >= 3*ny) & ...
	( peak(:,4) > (min(x)+0.2*dx ) & ( peak(:,4) < (min(x)+0.8*dx))));

peakint = peak(index,5);
peakpos = peak(index,1);	% test symetry around peaks and between them
lp = length(index);
peakpos = [ peakpos ; round(peakpos(1:(lp-1))+diff(peakpos)/2) ];
peakint = [ peakint ; round(peakint(1:(lp-1))+diff(peakint)/2) ];
[ peakpos, i ] = sort(peakpos);
peakint = peakint(i);
cencrit = ones(1,lp);

% now extend y signal

ny = [ ones(ly,1)*y(1) ; y ; ones(ly,1)*y(ly) ];

wy = min(max(peak(:,2)),max(peak(:,3)));	% should look around peakpos within wy
nx = [];
for i = 1:length(peakpos)		% working on indexes
	tx = (max(1,peakpos(i)-ceil(wy/2))):(min(ly,peakpos(i)+ceil(wy/2)));
	if isempty(nx)
		nx = tx(:);
	else
		nx = [ nx ; tx(:) ];
	end
end
if length(peakpos) < 2 | isempty(nx)
	disp('Wait during symetry serach...')
	nx = find( x > (min(x)+0.2*dx) & x < (min(x)+0.8*dx));
	cencrit = 0*nx;
end
t=nx;
fprintf(1,'Search symetry : ');
for i = 1:length(nx)
	xo = nx(i);
	lenx = ceil(max(xo,ly-xo)/2);
	xax = (xo-lenx):(xo+lenx);
	ny2 = ny-min(ny(xax+ly));
	ny2 = ny2/sum(ny2(xax+ly));
	cencrit(i) = sum( (ny2(xax+ly) - ny2(2*xo-xax+ly)).^2 );
end

% now get minimum of cencrit, in large 20%-80% zone
% plot(nx,t,nx,cencrit)
cencrit = cencrit/mean(cencrit);
[cencriteria, centre] = min(cencrit);
centre = x(nx(centre));
centresav = centre;

% look around gravity centre

gravit = sum(y.*x)/sum(y);
dx = abs(centre-gravit);
i = find( x(nx) > centre-3*dx & x(nx) < centre+3*dx);

if length(i) > 4
	nx = nx(i);
	cencrit2 = cencrit(i);
	[cencriteria2, centre] = min(cencrit2);
	dc = diff(cencrit2);
	centre = x(nx(centre));
	cencriteria = cencriteria2/max(abs(dc));
end

if isempty(centre)
	if ~isempty(centresav)
		centre = centresav;
	else
		disp('Gravity centre was used')
		centre = gravit;
	end
end
fprintf(1,'%f \n',centre);

