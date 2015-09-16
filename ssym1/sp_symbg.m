function [newy, newerr, centre,p,rightoverleft,z,correction] = sp_symbg(x,y,err,peakanalysis, centre, peakmat,analysispart,order);
% [newy, newerr, centreindex, p] = sp_symbg(x, y, ,err, peakanalysis, centre, peakmat,analysispart,order)
% This function makes signal background symetric in intensity
% p and centre are fit polynome and signal centre.
% (peakanalysis) is the signal analysis index vector, beeing :
%      3 for 'background' part.
%      1 for 'shutter/apparatus' zones.
%      2 for significative peaks.
% (analysispart), with same signification, is the part of signal to make symetric (default is background).
% Warning : underlying slopes/broad peaks are affected.
%           this operation must be done only if your signal is intrinsically symetric

% E.Farhi 04/98 (manuf@ldv.univ-montp2.fr)
% uses : mf_flsqr, sp_gappf, quadrat

% get signal between boundaries
% search centre
% compute AS/S
% fit AS/S (2nd order)
% modify y for symetry

% first extract data not in shutter/apparatus zones and not in significant peaks
% [ index left_width right_width  max_pos Intensity Width ]

x=x(:); y=y(:);
if ~exist('order') | isempty(order)
	order = 2;
end
if ~exist('analysispart') | isempty(analysispart)
	analysispart = 0;
end
ly = length(y);

newy = y; p = []; newerr = err;
if ~exist('peakmat') | ~exist('peakanalysis') | isempty(peakmat) | isempty(peakanalysis)
	[shutterzones, peak, peakanalysis, levely, sp_type, pmin,pmax] = sp_gappf(x,y);
end

if ~exist('centre') | isempty(centre)
	centre = sp_cent(x,y,peak);;
end

% get centre in indexes

lix = find(x <= centre);
rix = find(x >= centre);
centresav = centre;
if isempty(lix) | isempty(lix)
	disp('warn : center not found in X axis !!')
	return
end
lix = lix(length(lix));
rix = rix(1);
[maxx,i] = min([ max(abs(x(rix) - centre)) max(abs(x(lix) - centre)) ]);
if i == 1
	centre = rix;
else
	centre = lix;
end


% now get half signal

index_to_compare = min(centre,ly-centre)-1;
index_to_compare = 1:index_to_compare;
bip = find(y(centre-abs(-index_to_compare)) == 0);

zx = index_to_compare;
if ~isempty(bip)
	index_to_compare(bip) = []; % suppress zero data before division
end
rightoverleft = y(centre+index_to_compare) ./ y(centre-abs(-index_to_compare));

di = diff(rightoverleft); di = [ di(1) ; di(:) ];
rapportok = (abs(di) < median(abs(di)) & ~isnan(rightoverleft) & ~isinf(rightoverleft) & rightoverleft > .25 & rightoverleft < 4);

if length(analysispart) == 1 & analysispart
	toanal = (peakanalysis(centre+index_to_compare) == analysispart | peakanalysis(centre-abs(-index_to_compare)) == analysispart);
else
	toanal = rapportok;
end

toanal = find(toanal & rapportok);

% theoretical result should be 1 on background.
if length(toanal) == 0
	disp('warn : no points for background symetrization')
	p = [];
	return
end
p = polyfit(toanal,rightoverleft(toanal),order);
dp = ones(size(p));
p(length(p)) = 1;
dp(length(p)) = 0;
p = mf_flsqr(toanal,rightoverleft(toanal),[],p,dp,'polynomial',[]);
% fit with no output
% make division equals 1 on symetry centre (a=1)

correction = polyval(p,1:centre);
y(centre:-1:1) = y(centre:-1:1).* correction(:);
err(centre:-1:1) = err(centre:-1:1).* correction(:);
newy = y; newerr = err;
