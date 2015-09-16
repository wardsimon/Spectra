function s=rebin(s, dx, method)
% function s=rebin(s, dx)
%
% SPED1D/REBIN Rebins a spec1d spectrum.
%
% dx = scalar: specifies bin widths.
% dx = vector: specifies the bin centers.
% dx = sped1d: uses x-values of dx as bin centers
%
% The width of boundary bins is twice the distance to neighbour bin.
%
% method='average' : New x values are average within bins
%                        (weighted by error)
%        'force'   : New x values are as specified in dx
%        'interp'  : New x-values as specified by dx, y and e
%                        interpolated from 'average'-x, 
%
% HMR 20.11.2000

if nargin<1
  error('No input to rebin')
end

%----- Get bin centers from dx
nps=length(s.x);
if nargin<2
%  error('Need interface to select bins')
  dx=1*mean(s.x(2:end)-s.x(1:end-1));
  dxi=input(['Input bin width [' num2str(dx) '] ']);
  if isreal(dxi)
    dx=dxi;
  end
end

if isa(dx,'spec1d')
  xbin=getfield(dx,'x');
elseif length(dx)==1
  xbin=min(s.x):dx:max(s.x);
elseif isa(dx,'char') & strcmp(dx,'auto')
  dx=2*mean(s.x(2:end)-s.x(1:end-1));  
  xbin=min(s.x):dx:max(s.x);
else
  xbin=dx;
end

% Establish the bin boundaries
xbin=xbin(:);
xbnd=[xbin(1)-(xbin(2)-xbin(1))/2
     (xbin(1:end-1)+xbin(2:end))/2
      xbin(end)+(xbin(end)-xbin(end-1))/2];
% Fill each bin in turn
x=[];
y=[];
e=[];
for i=1:length(xbnd)-1
  j=find((s.x>=xbnd(i)) & (s.x<xbnd(i+1)));
  if ~isempty(j)
     % rebin uses simple addition of points.
     % One should be carefull if counts have different monitors.
    x=[x;mean(s.x(j))];
%    x=[x;sum(s.x(j)./s.e(j))/sum(1./s.e(j))];
    y=[y;sum(s.y(j))];
    e=[e;sqrt(sum(s.e(j).^2))];
  end
  
end

s.x=x;
s.y=y;
s.e=e;


   
if strcmp(method,'forced')
   s.x=xbin;
   s.y=y;
   s.e=e;
elseif strcmp(method,'interp')
   s=interpolate(s,xbin);
end