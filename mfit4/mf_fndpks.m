function [peak,pmax,sy,ny,sdy,ndy,sddy,nddy,pmin] = mf_fndpks(x,ys,thr, pos,width,pmin,pmax)
%Quick search of peaks in a signal.
%Syntax: [peakmat,pmax,sy,ny,sdy,ndy,sddy,nddy,pmin] = mf_fndpks(x,y,{threshold=1, pos,width,,pmin,pmax}) or mf_fndpks(y)
%
% This function searches peaks into an (x,y) signal, using a specified
% 'threshold' in noise units, usually around 1 ; few peaks for a high value.
% 'peakmat' is an array of rows for each maximum identified
%          [ index left_width right_width  max_pos Intensity Width MinIdxBefore MinIdxAfter ].
% 'pmax' is a Dirac type peaks signal.
% Averaged signals (including derivatives) and noises are also returned.
% opt : pos, width = part of 'y' to be scanned for peaks around pos.
%       if pmin and pmax are given, a faster guess is done.
%       smoothing order is 2.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: quick search of peaks in a signal.

% Part of 'Spectral tools'. E.Farhi. 07/96
% uses : diff.m
%  mf_smooth.m


if (nargin < 1)
  error('usage : [peaks,pmax,sy,ny,sdy,ndy,sddy,nddy] = mf_fndpks(x,y,{threshold=1}) or peaks(y)');
end

if nargin < 7, pmax = ''; end
if nargin < 6, pmin = ''; end
if nargin < 5, width = ''; end
if nargin < 4, pos = ''; end
if nargin < 3, thr = ''; end
if nargin < 2
  ys=x;
  x=1:length(x);
end

if isempty(thr)
  thr = 1;
end

if isempty(pos) | isempty(width)
  nd = length(x);
  pos = x(ceil(nd/2));
  width = abs(x(nd)-x(1));
end

logflag = 0;

if thr == 0, thr = 0.0001; end

y=ys(:);

nd = length(y);
dx = x(2) - x(1);

%if any(diff(x) ~= dx) % step is not constant, need to interpolate
% nx = linspace(x(1), x(nd), round((x(nd)-x(1))/dx));
% y = interp1(x,y,nx);
% x=nx;
% plot(x,y)
%end

t=find( (x>=(pos-width)) & (x<=(pos+width)) );    % index zone around pos
t = t(:);
y = y(t); y=y(:);
x = x(t); x=x(:);
nd = length(y);
miny = min(y);
if miny == 0
  y=y+min(abs(y(find(y~=0)))/10);
end

if (max(y)/median(y) > 1e2)
  logflag = 1;
  y = log(y-min(y)+1);
end

if isempty(pmin) | isempty(pmax)
  sy = mf_smooth (y);
  dy = diff (sy); dy(nd) = dy(nd-1);
  sdy = mf_smooth (dy);
  ddy = diff (sdy); ddy(nd) = ddy(nd-1);
  sddy = mf_smooth (ddy);

  t=(abs(y) <= median(abs(y)));   % noises computation
  t = t(:);
  nt = sum(t);
  ny = std(t.*y - t.*sy);
  t=(abs(dy) <= median(abs(dy)));
  ndy = std(t.*dy - t.*sdy);
  t=(abs(ddy) <= median(abs(ddy)));
  nddy = std(t.*ddy - t.*sddy);

  t = dy(1:(nd-1)).*dy(2:nd);   % for 1st derivative sign changes.
  t(nd) = sdy(nd)^2;
  t = t(:);

  pmax = (abs(dy) <= ndy/thr) | (t <= (ndy/thr)^2); % extrema

  pmin = ( pmax & (sddy > nddy*thr));
  pmax = ( pmax & (sddy < -nddy*thr));
end

% now alterning maxima and minima ...

pmaxsav = pmax;
peaklist=find(pmax);
index =1;   % index in original peak list (from pmax)
ip = 1;   % index in final peak table
peak = [ 1 1 1 0 0 0 0 0 ];
pmax = 0*y;

while (index <= length(peaklist))  % scanning all peaks in pmax
  peakindex = peaklist(index);
  t = find(pmin(1:peakindex));  % search min before
  if (~isempty(t))
    minbefore = t(length(t));
  else
    minbefore = 1;
  end
  t = find(pmin(peakindex:nd)); % search min after
  if (~isempty(t))
    minafter = t(1) + peakindex -1;
  else
    minafter = nd;
  end
  t = minbefore:minafter;
  [dummy,t] = max(y(t));  % find real max beween 2 min.
  if (isempty(t))
    index = index+1;      % no max found
    % trying next in pmax
  else
  max_pos = t(1) + minbefore -1;      % max found
  t = min(y(minbefore),y(minafter));
  if ((y(max_pos) - t) <= ny) % this is a noise peak.
% fprintf(1,'[%i] %i %i %i : %i %i \n',max_pos,y(max_pos),t,ny,y(minbefore),y(minafter));
    if (y(minbefore) > y(minafter)) % equivalent minima
      pmin(minbefore) = 0;
    else      % higher minimum removed
      pmin(minafter) = 0;
    end
    if (minafter == nd)
      break;      % finish
    else
      if (minbefore == 1)
        index = index +1;
      end  % trying again

    end
  else        % real peak
    if (~isempty(find((peak(:,1) - max_pos) == 0)))
      index = index + 1;    % already found
    else
    peak(ip,1) = max_pos;
    t = min(y(minbefore),y(minafter));
    t = 0.1*y(max_pos) + 0.9*t; % 10 percent hight for peak
    halfafter = find(y(max_pos:minafter) < t);
    if (isempty(halfafter))
      halfafter = minafter;
    else
      halfafter = halfafter(1) + max_pos -1;
    end
    halfbefore = find(y(minbefore:max_pos) < t);
    if (isempty(halfbefore))
      halfbefore = minbefore;
    else
      halfbefore = halfbefore(length(halfbefore)) + minbefore -1 ;
    end

    peak(ip, 2) = max(1,max_pos - halfbefore);
    peak(ip, 3) = max(1,halfafter - max_pos);
    t = halfbefore:halfafter;
%disp([ max_pos length(t) ])
% might be done with polyfit but this is better...
%   [ smoothed, p] = mf_interp(x(t),y(t),x(max_pos),max(ceil(length(t)/2),3),2);
    p = polyfit(x(t),y(t),2);

    S = x(max_pos);
    W = abs(x(halfafter) - x(halfbefore));
    In = y(max_pos);

    if (length(p) == 3)
      a=p(1); b=p(2); c=p(3); % 2nd order interpolation
      delta = (b*b - 4*a*c);
      if ((a < 0) & (delta >0))
          S=-b/2/a;
        W=abs(sqrt(delta)/a);
        In=a*S*S+b*S+c;
      else
        t = minbefore:minafter;
%       [ smoothed, p] = mf_interp(x(t),y(t),x(max_pos),max(ceil(length(t)/2),3),2);
        p = polyfit(x(t),y(t),2);

        if (length(p) == 3)
          a=p(1); b=p(2); c=p(3); % 2nd order interpolation
          delta = (b*b - 4*a*c);
delta = -1;
          if ((a < 0) & (delta >0))
            S=-b/2/a;
            W=abs(sqrt(delta)/a);
            In=a*S*S+b*S+c;
          else
            p = [];
          end
        end
      end
    else
      p = [];
    end
    if S > max(x) | S < min(x) | W > abs(max(x)-min(x))
      p = [];
    end
    if isempty(p)     % interpolation not performed
      S = x(max_pos);
      W = abs(x(halfafter) - x(halfbefore));
      In = y(max_pos);
    end

    t = [ minbefore minafter ];
    tmin = t;
    wx = find(x > S+W);
    if ~isempty(wx), t = [ t wx(1)]; end
    wx = find(x >= S-W);
    if ~isempty(wx), t = [ t wx(1) ]; end
    t = abs(t - max_pos) * 3;
    t = max(t);
    t = (max_pos -t):(max_pos+t);
    t = t(find(t>=1 & t<= nd));
    t = union(t,tmin);
    t = min(y(t)); % compared with nearer min peak
    peak(ip,4) = S;
    peak(ip,5) = In - t;
    peak(ip,6) = W;
    peak(ip,7) = minbefore;
    peak(ip,8) = minafter;
    ip = ip + 1;
    index = index + 1;
    pmax(max_pos) = 1;
    end % from if (~isempty(find((peak(:,1) - max_pos) == 0)))
  end
  end
end

if (logflag)
  peak(:,5) = exp(peak(:,5)) + miny -1;
end
