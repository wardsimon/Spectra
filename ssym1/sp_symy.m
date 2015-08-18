function [newx,newy,newerr,tosym] = sp_symy(x,y,err,centre,peak,peakanalysis,option)
% [newx,newy,newerr] = sp_symy(x,y,err,centre,peak,peakanalysis,option)
% used after sp_symbg (optional) and/or sp_symx.
% options : boolean vector for options (default is [ 1 1 1 2 ])
%        * symetrization for background (1=yes).
%        * symetrization for apparatus (1=yes).
%        * symetrization for signal (1=yes).
%        * order

% get spectra centre
% get left and right part
% make mean x axis -> newx (interpolation), then interpolate y on it
% compute newy locally or globally.

% E.Farhi 04/98 (manuf@ldv.univ-montp2.fr)
% uses : sp_symbg

if ~exist('option') | isempty(option)
	option = [ 1 1 1 2 ];
end

lix = find(x <= centre);
rix = find(x >= centre);

maxx = min(max(abs(x(rix) - centre)), max(abs(x(lix) - centre)));
newx = linspace(-maxx, +maxx, length(y));
disp('Working...')
tosym = 0*y;

order = option(4); newerr = [];

if option(2)
	tosym(find(peakanalysis == 1)) = 1;
end
if option(3)
	tosym(find(peakanalysis == 2)) = 1;
end
if option(1) | all(option == 0)
	tosym(find(peakanalysis == 3)) = 1;
	if option(1), analysispart = 3; else analysispart = 0; end
	[y, newerr, dummy1, dummy2] = sp_symbg(x, y, err, peakanalysis, centre, peak,analysispart,order);
	if ~isempty(dummy2)
		disp('Background symetrization performed')
	end
end
if option(2) == 0 & option(3) == 0
	newy = y;
	newx = x;
	newerr = err;
	tosym = peakanalysis;
	return
end

newy = interp1(x,y,newx+centre,'spline');
if isempty(newerr)
	newerr = interp1(x,err,newx+centre,'spline');
end
tosym = interp1(x,tosym,newx,'nearest');

lix = find(newx < 0); lix = lix(:);
rix = find(newx > 0); rix = rix(:);	% should be of same lenghts
flix = flipud(lix);	% revert : start from center

ny = (newy(flix) + newy(rix))/2;	% mean
newyo = newy;
tl = tosym(flix);
tr = tosym(rix);

toset = find(tr == 1 & tl == 1);
toset = toset(:);

newy(flix(toset)) =  ny(toset);
newy(rix(toset)) = ny(toset);

apny = find(tosym == 1);
apy = find(peakanalysis == 1);
i = find(newy(apny) < min(y(apy)));

if ~isempty(i)
	newy(i) = min(y(apy));
end

tosym = interp1(x,peakanalysis,newx,'nearest');
disp('Full symetrization performed')
