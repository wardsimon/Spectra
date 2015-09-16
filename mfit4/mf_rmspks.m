function [ny,npeak] = mf_rmspks(x,y,peak,levely)
% [ny,npeak] = mf_rmspks(x,y,peak) : remove spikes

toremove = []; npeak = [];
nx = x; ny = y;

if nargin < 4, levely=[]; end
if nargin < 3, peak=[]; end
nny = [];

if isempty(peak)
  [peak,pmax,sy,nny,sdy,ndy,sddy,nddy,pmin] = mf_fndpks(x,y,0);
end

if isempty(levely)
  if ~isempty(nny)
    levely = 4*nny;
  else
    levely = min(y);
  end
end

if (max(peak(:,2) .* peak(:,3)) > 2)
  i = find( peak(:,2) .* peak(:,3) <= 2);
  toremove = [toremove ; peak(i,:) ];
  peak(i,:) = []; % remove spikes

  i = find( (y(peak(:,1)) < levely) | (peak(:,5) < levely) );
  toremove = [toremove ; peak(i,:) ];
  peak(i,:) = []; % remove small peaks

  i = find( min(peak(:,2), peak(:,3)) == 1);
  toremove = [toremove ; peak(i,:) ];
  peak(i,:) = []; % remove 1 width peaks
end

npeak = peak;

for ip = 1:size(toremove,1)
  ny(toremove(ip,7):toremove(ip,8)) = min(y);
end
