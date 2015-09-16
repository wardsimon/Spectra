function [X,Y,Z,h]=mesh3(x,y,z,method,dims)
% mesh3 : plot a 3d mesh from (x,y,z) points (see plot3)
% [X,Y,Z,h] = mesh(x,y,z,{method})
% plots a 3D surface from individual (x,y,z) points.
% method is a string that can contain options :
% * griddata method
%    cubic, linear, nearest, v4
% * graph3d option
%    plot3, mesh, surf, surfc, surfl, fill3,
%    contour, contour3, contourf, stem3
% * shading option and others
%    flat, interp, faceted
%    clabel (for contour's only)
%    labels
% ex: mesh3(x,y,z,'cubic+surf+interp');
% default is linear+mesh

if (nargin == 1)
  error('mesh3 require 2 or 3 parameters (vectors)')
end
if (nargin == 2)
  if (ischar(y))
    error('mesh3 require at least x and z vectors');
  end
  z = y;
  y = x;
end

if (nargin == 3) & (ischar(z))
  method = z;
  z = y;
  y = x;
end

if nargin < 3, dims=[]; end
if isempty(dims)
    dims = [33 33];
end

if nargin < 4,
  method = '';
end
if isempty(method)
  method = 'linear+mesh';
end

x = x(:); y=y(:); z=z(:);
xlen=prod(size(x));
if (length(y) ~= xlen) | (length(z) ~= xlen)
  error('Vectors x,y and z should have the same length');
end

% sorting data along x
[x,xsorti] = sort(x);
y = y(xsorti);
z = z(xsorti);
clear xsorti

% resampling of data for mesh
xmin=min(x);
ymin=min(y);
xmax=max(x);
ymax=max(y);

xlin=linspace(xmin,xmax,dims(1));
ylin=linspace(ymin,ymax,dims(2));

[X,Y] = meshgrid(xlin,ylin);
Z = [];
clear xlin ylin

if (findstr(method,'v4'))
  gridmethod = 'v4';
elseif (findstr(method,'cubic'))
  gridmethod = 'cubic';
elseif (findstr(method,'nearest'))
  gridmethod = 'nearest';
else
  gridmethod = 'linear';
end

if isempty(findstr(method,'plot3'))
  Z = griddata(x,y,z,X,Y,gridmethod);
  clear x y z
  method = [ method '    ' ];
  C = [];
  if (findstr(method,'contour3'))
    [C,h] =contour3(X,Y,Z);
  elseif (findstr(method,'contourf'))
    [C,h]=contourf(X,Y,Z);
  elseif (findstr(method,'contour'))
    [C,h]=contour(X,Y,Z);
  elseif (findstr(method,'surfc'))
    h=surfc(X,Y,Z);
  elseif (findstr(method,'surfl'))
    h=surfl(X,Y,Z);
  elseif (findstr(method,'surf'))
    h=surf(X,Y,Z);
  elseif (findstr(method,'stem3'))
    h=stem3(X,Y,Z);
  else
    h=mesh(X,Y,Z);
  end

  if (findstr(method,'flat'))
    shading flat
  elseif (findstr(method,'interp'))
    shading interp
  elseif (findstr(method,'faceted'))
    shading faceted
  end

  if ~isempty(C) & findstr(method,'clabel')
    clabel(C,h);
  end

  if (findstr(method,'labels'))
    xlabel(inputname(1));
    ylabel(inputname(2));
    zlabel(inputname(3));
    title([ inputname(3) ':' method ]);
  end

else
  h=plot3(x,y,z,'.');
end

axis tight
%shading interp
rotate3d on
if (nargout == 0) X=method; end

