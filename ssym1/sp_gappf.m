function [shutterzones, peak, peakanalysis, levely, sp_type, pmin,pmax,pass1,pass2,peak1] = sp_gappf(x,y,levely,pmin, pmax,selected)
% [shutterzones, peak, peakanalysis, levely, sp_type, pmin,pmax] 
%     = sp_gappf(x,y,{levely,pmin,pmax})
% This function looks for the instrumental apparatus function into a sprectra.
% Input : (x,y) signal and optionnaly detection level and peak vectors
% Output: shutterzones  a [start end] matrix for shutter zones.
%         shutter       1 when shutter is closed, 0 for background, 2 = signal.
%         centre        probable spectrum centre (if any).
%         peak          list of peaks [ index left_width right_width  max_pos Intensity Width MinIdxBefore MinIdxAfter ]
%         levely        detection level used for signal analysis.
%         sp_type       spectrum type 'flips' or 'peaks'.
%         cencriteria   symetry for < 1e-1

% E.Farhi 04/98 (manuf@ldv.univ-montp2.fr)
% uses : mf_fndpks

flagout = 1;

if ~exist('levely')
	levely = [];
end

if nargin > 3 & (~exist('selected') | isempty(selected))
	flagout = 0;
end

if ~exist('pmin'), pmin = []; end
if ~exist('pmax'), pmax = []; end

i = find(y <= 0);
if ~isempty(i)
	y(i) = min(y(find(y>0))/10);
end
x=x(:); y=y(:);
xs= x; ys = y;
ly = length(y);

if ~exist('selected') | isempty(selected)
	selected = ones(size(y));
end

shutterzones = [];
sp_type = 'flips'; peakanalysis = [];

% detect a kind of noise threshold... ----------------------------

levely0 = median(y);
%fprintf(1,'Signal base threshold is : %f\n',levely0);

medmat = [];
my = (median(y) + mean(y))*2;
n = 100;
while (n > 0)
	myo = my;
	my = (median(y) + mean(y));
	medmat = [ medmat my ];
	i = find(y < my);
	if isempty(i) & flagout 
		disp('No more points');
	end
	y = y(i);
	n = n - 1;
	if ((myo /my) - 1) < 1e-3,
		break
	end
end

levely1 = medmat(length(medmat));

%fprintf(1,'Signal first threshold is : %f\n',levely1);
histy1 = y(find(y<=median(y)));
if isempty(histy1)
	histy1 = y(find(y<=mean(y)));
end
[histy1,histx1] = hist(histy1);
[dummy, levely2] = min(histy1); % find minimum of y distribution under threshold
levely2 = levely2(1);
if (levely2 == 1) | (levely2 == 10)
	if flagout
		disp('There seems to be no shutter zones...');
	end
	levely2 = levely1;
else
	p = polyfit(histx1((levely2-1):(levely2+1)), histy1((levely2-1):(levely2+1)),2);
	a=p(1); b=p(2); c=p(3);	% 2nd order interpolation

	if (a ~= 0)
		levely2 =-b/2/a;
	else
		if flagout 
			disp('Can''t find minimum of distribution');
		end
		levely2 = histx1(levely2+1);
	end 
end
%fprintf(1,'Signal 2nd threshold is : %f\n',levely2);

y = ys;
x = xs;

if isempty(levely)
	levely = min([ levely0 levely1 levely2 ]); 	% shutter zones are under this level
end

if flagout 
	fprintf(1,'Low level is %.2g in signal\n',levely);
end

% searching peaks in spectrum ----------------------------

[peak,pmax,sy,ny,sdy,ndy,sddy,nddy,pmin] = mf_fndpks(xs,y,0,[],[],pmin, pmax); % max number of peaks
% [ index left_width right_width  max_pos Intensity Width ]

peakall = peak;

% clearing 'bad' peaks ----------------------------

[ny, peak] = mf_rmspks(x,y,peak,levely);

peaknobad = peak;

% analysing shutter zones by asymteric criterium ------------

peak = peakall;
wasym = log(peak(:,2)./peak(:,3));
w = min(peak(:,2),peak(:,3))*2;
Iasym = log (y(max(peak(:,1) - w,1)) ./ y(min(peak(:,1) + w,ly)));
% Iasym and wasym are <0 for left shutter edge, >0 for right.

dummy = find(y(peak(:,1)) > levely & selected(peak(:,1)));	% remove small peaks
peak  = peak(dummy,:);
wasym = wasym(dummy,:);
Iasym = Iasym(dummy,:);

dummy = find(wasym.*Iasym > 0);		% get asymteric peaks
peak  = peak(dummy,:);
wasym = wasym(dummy,:);
Iasym = Iasym(dummy,:);

[t1, dummy] = sort(wasym.*Iasym);	% sort by descending asymtery
dummy = flipud(dummy);
peak  = peak(dummy,:);
wasym = wasym(dummy,:);
Iasym = Iasym(dummy,:);

shutinterval = mean(diff(sort(peak(:,1))));
t1 = [ peak(:,1), wasym Iasym wasym.*Iasym y(peak(:,1)) sign(wasym) ];

pass1 = []; pass2 = []; peak1 = peak;

i = 1;
while i <= size(t1,1)	% look for asymetric edges
	pos1 = t1(i,1);
	typ1 = t1(i,6);
	int1 = t1(i,5);
	notfound = 0;
	if typ1 == 1	% other edge should be : opposite sign, quite near
		i2 = find(t1(:,1) > pos1 & t1(:,6) == -1 & t1(:,1) - pos1 < 3*shutinterval);
	else
		i2 = find(t1(:,1) < pos1 & t1(:,6) ==  1 & pos1 - t1(:,1) < 3*shutinterval);
	end

	if ~isempty(i2)
		t2 = t1(i2,:);
		[dummy, i2] = sort(abs(t2(:,1) - pos1)); % search nearest peak
% interval should contain in each interval at leat two minima and one peak
		geti = 1;
% fprintf(1,'\n[ %i ] ',pos1);
		while geti <= min(length(i2),1)
			pos2 = t2(i2(geti),1);
			nmin = find(pmin);
			nmin = find(nmin > min(pos1,pos2) & nmin < max(pos1,pos2));
			nmax = find(peaknobad(:,1) >= min(pos1,pos2) & peaknobad(:,1) <= max(pos1,pos2)); 
% fprintf(1,': %i %i %i ',pos2,length(nmin),length(nmax));
			if length(nmin) >= 2 & length(nmax) >= 1
				break
			else
				geti = geti+1;
			end
		end
		if geti <= min(length(i2),1)	% previous condition is ok
			i2 = i2(geti);
			pos2 = t2(i2,1);
			typ2 = t2(i2,6);
			int2 = t2(i2,5);

			toadd = sort([ pos1 pos2 ]);
			pos1 = min(toadd); pos2=max(toadd);
			if min(y(pos1:pos2)) < 2*levely
% is it really a shut ?
				if isempty(pass1)
					pass1 = toadd;
				else
					pass1 = [ pass1 ; toadd ];
				end
			elseif flagout
				fprintf(1,'Not min enough in [%i:%i]',pos1,pos2);
				notfound = 1;
			end
		elseif flagout
			fprintf(1,'No min/max for [%i]',pos1);
			notfound = 1;
		end
	elseif flagout
		fprintf(1,'No edge for [%i]',pos1);
		notfound = 1;
	end
	if notfound
		if peak(i,2) > peak(i,3)
			pos2 = min(ly,pos1 + 3*peak(i,3));
		else
			pos2 = max(1,pos1 - 3*peak(i,2));
		end
		toadd = sort([ pos1 pos2 ]);
			pass1 = [ pass1 ; toadd ];
		fprintf(1,' gather with [%i]\n',pos2);
	end
	i = i + 1;
end

if ~isempty(pass1)

[dummy, i] = sort(pass1(:,1));
pass1 = pass1(i,:);

shuti = 1;	% shutter window index in pass1
while shuti <= size(pass1,1) 	% group shutter windows and intersect.
	b1 = pass1(shuti,1);
	b2 = pass1(shuti,2);
	shuttable = [ b1 +1 ; b2 -1 ];
	shuti = shuti + 1;
	groupi = shuti;
	while groupi <= size(pass1,1)
		b1 = pass1(groupi,1);
		b2 = pass1(groupi,2);
		if ~isempty(find(shuttable == b1)) | ~isempty(find(shuttable == b2))
			shuttable = [ shuttable ; b1 +1 ; b2 -1 ];
			groupi = groupi +1;

		else
			groupi = size(pass1,1)+1;
		end
	end
	[dummy, sorti] = sort(shuttable(:,1));
	shuttable = shuttable(sorti,:);

% each possible shutter zone is now isolated and sorted.
% get intersection
	b1 = max(shuttable(find(shuttable(:,2) == 1),1));
	b2 = min(shuttable(find(shuttable(:,2) == -1),1));
	if abs(log(levely/max(y(b1:b2)))) > 2
		if y(b1) > levely
			b1 = b1-abs(b1-b2);
		end
		if y(b2) > levely
			b2 = b2+abs(b1-b2);
		end
		if isempty(pass2)
			pass2 = [ b1 b2 ];
		else
			pass2 = [ pass2; b1 b2 ];
		end
	end
end

end % if ~isempty(pass1)

if ~isempty(pass2)

[dummy, i] = unique(pass2(:,1));
pass2 = pass2(i,:);

% now we remove those boundaries from peaknobad matrix

peak = peaknobad;
for pki = reshape(pass2,1,prod(size(pass2)))
	i = find(peak(:,1) == pki);
	if ~isempty(i)
		peak(i,:) = [];
	end
end
peaknobad = peak;

% then we modify shutterzones with edges width

peak = peakall;
lp2 = size(pass2,1);
for shuti = 1:lp2
	b1 = max(1,pass2(shuti,1));
	b2 = min(length(y),pass2(shuti,2));
	if (shuti == lp2) | (shuti < lp2 & (pass2(shuti+1,1) > b2))
		[my, bm] = max(y(b1:b2));
		[my, c1] = min(y(b1:(b1+bm-1))); c1 = c1 + b1 -1;
		[my, c2] = min(y((b1+bm-1):b2)); c2 = c2 + b1 + bm -2;
		b1 = max(1,c1);
		b2 = min(ly,c2);

		if isempty(shutterzones)
			shutterzones = [ b1  b2 ];
		else
			shutterzones = [ shutterzones ; b1  b2 ];
		end
	end
end

end % if ~isempty(pass2)

% now look for peak with minimum width :  apparatus function
% also look if there are multiple thin peaks (for SHR for instance)

peak = peaknobad;
lp = size(peak,1);

nb_shut = size(shutterzones,1);
if (nb_shut == 0) & (lp >= 1) | (levely2 == levely1)
	sp_type = 'peaks';
	if flagout
		disp('Trying peak discrimination...');
	end
	sy = mf_smooth(y);
	ny = std(y-sy);

	i = find((y(peak(:,1)) < levely + 5*ny) & (peak(:,5) < 5*ny));
	peak(i,:) = [];
	if (length(peak(:,5)) < 1)
		if flagout
			disp('Not enough peaks in signal.');
		end
		shutterzones = [];
		shutter = 0*y;
	else
% [ index left_width right_width  max_pos Intensity Width ]
%	[ peak(:,[1 4 5 6]) y(peak(:,1)) ]
	[hy,hx] = hist(peak(:,5),ceil(size(peak,1)));
	[a,b] = min(hy);
	hx = hx(b); % first minimum level in peak intensities
	ipeak = find(peak(:,5) > hx); % & peak(:,4)-2*peak(:,6) > min(x) & peak(:,4)+2*peak(:,6) < max(x)); 	% can add | y(peak(:,1)) > levely+hx)

	if ~isempty(ipeak)
		mpeak = peak(ipeak,1);
		[wmpeak, sorti] = sort(peak(ipeak,6));
		xmpeak = x(mpeak); xmpeak = xmpeak(:);
		ipeak = wmpeak+xmpeak;
		i = find(ipeak > ly); ipeak(i) = ly;
		wmpeak = round(mf_interp(x,1:ly,ipeak)); 

		wmpeak = wmpeak(:) - mpeak;
		mpeak = mpeak(sorti);
		lastw = wmpeak(1);
% [ mpeak wmpeak ]
		for i = 1:length(wmpeak)
			max_pos = mpeak(i);
			if abs(log(wmpeak(i)/lastw)) < .5
				shutterzones = [ shutterzones ; max(1,(max_pos-wmpeak(i))) min(ly,(max_pos+wmpeak(i))) ];
				lastw = wmpeak(i);
			end
		end
	end
	end % if (length(peak(:,5)) < 3)
end
nb_shut = size(shutterzones,1);
if flagout
	fprintf(1,'Found %i apparatus zones\n',nb_shut);
end

shutter = 0*y;
for i=1:size(shutterzones)
	b1 = shutterzones(i,1);
	b2 = shutterzones(i,2);
	shutter(b1:b2) = 1;
end

peak = peaknobad;	% add eventual forgotten apparatus peaks
i = find( (shutter(peakall(:,1)) == 1));
if ~isempty(i)
	peak = [ peak ; peakall(i,:) ];
	[j,k] = unique(peak(:,1));
	peak = peak(k,:);
end
peaknobad = peak;

% separate data into background, apparatus and peak parts.

np = size(peak,1); % number of significant peaks
peakanalysis = 3*ones(size(x));
for i = 1:np
	peakanalysis(max(1,peak(i,1)-2*peak(i,2)):min(ly,peak(i,1)+2*peak(i,3))) = 2;
end
peakanalysis(find(shutter)) = 1;

peak = peakall;
