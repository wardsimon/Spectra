function [xi,yi,zi]=mapplotmodifiziert(s,extvar,smoothe,xi,yi,scale,laenge)
%
% function [xi,yi,zi]=mapplot(s,extvar,smoothe,xi,yi,scale)
%
% SPEC1D/mapplot creates 2D map of an array of spectra
%
% extvar		x-variable that separates the spectra
% smoothe	each spectrum can be convluted by a gaussian
% xi			x-bins for map (leave [] for automatic choice)
% yi			y-bins for map (leave [] for automatic choice)
% scale		if scale='log'
%
% if no output is requested, a colorplot of the data is made.
%
% This function is a visualisation tool - not to be used for data manipulation.
%
% HMR 20.11.2000

if nargin<2
   extvar=1:length(s);
end
if nargin<3
   smoothe=[];
end

x=[];
y=[];
z=[];
for n=1:length(s)
   if ~isempty(smoothe)
      s(n)=smooth(s(n),smoothe);
   end
   y=[y;getfield(s(n),'x')];
   z=[z;getfield(s(n),'y')];
   if size(extvar)==[length(getfield(s(n),'x')) length(s)]
      x=[x;extvar(:,n)];   
   else
      x=[x;extvar(n)*ones(size(getfield(s(n),'x')))];   
  end
end

if nargin<5
   xi=[];
   yi=[];
end
if nargin>5
   if strmatch('log',scale)
      z=log(z);
   end
end
if isempty(xi) | isempty(yi)
%   xi=min(x):(max(x)-min(x))/100:max(x);
%   yi=min(y):(max(y)-min(y))/100:max(y);
   xi=min(x):(max(x)-min(x))/200:max(x);
   yi=min(y):(max(y)-min(y))/200:max(y);
end
xi=xi(:)';
yi=yi(:);

[xi,yi,zi]=griddatamodifiziert(x,y,z,xi,yi,'cubic',laenge);

if nargout==0
  pcolor(xi,yi,zi)
  shading flat
  line(x,y,'linestyle','none','marker','.','color','k')  
end

return

if 1==0
   % Do the triangulation thing, which is complicated but might be 
   % pretier for more dispersieve modes. Not general and needs specific
   % recoding for you particular application.
   % HMR 7.12.2000
if isa(excit1,'spec1d')
   excit1=getfield(excit1,'y')
end

ss(1)=smooth(s(1),smoothe);
ynew=getfield(ss(1),'x');
znew=getfield(ss(1),'y');
xnew=field(1)*ones(size(ynew));
x=xnew;y=ynew;z=znew;
nprev=0;
ntot=length(xnew);
d=[];
for n=2:length(s)
   ss(n)=smooth(s(n),smoothe);
   xprev=xnew;
   yprev=ynew;
   zprev=znew;
   ynew=getfield(ss(n),'x');
   znew=getfield(ss(n),'y');
   xnew=field(n)*ones(size(getfield(ss(n),'x')));
   
if n==2   
   n1=max(find(yprev<excit1(n-1,1)));
   n2=max(find(ynew<excit1(n,1)));
   d1=delaunay([xprev(1:n1);xnew(1:n2)],[yprev(1:n1);ynew(1:n2)]);
   d1=d1+(d1<=n1)*nprev+(d1> n1)*(nprev+length(xprev)-n1);
   d2=delaunay([xprev(n1:end);xnew(n2+1:end)],[yprev(n1+1:end);ynew(n2:end)]);
   d2=d2+(d2<=(length(xprev)-n1))*(nprev+n1-1)...
        +(d2>(length(xprev)-n1))*(nprev+n1+n2-1);
else
   d1=delaunay([xprev;xnew],[yprev;ynew])+nprev;
   d2=[];
end
     
   y=[y;ynew];
   z=[z;znew];
   x=[x;xnew];
%   d=[d;delaunay(x(nprev+1:end),y(nprev+1:end))+nprev];
   d=[d;d1;d2];
   nprev=ntot;
   ntot=length(x);
end

clf
trisurf(d,x,y,z)
view(2)
%shading flat
shading interp
caxis([0 5])
axis([12 14.5 0 4])
end
