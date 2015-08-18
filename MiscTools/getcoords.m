function [nx,ny] = getcoords(handle)
% [nx,ny] = getcoords(handle) : show mouse position on current graph
% optional handle activates given figure
% clik on right button to end

if nargin > 0
  figure(handle);
end

button = -1;
nx = []; ny = [];
while button ~= 3
  [x,y,button] = ginput(1);
  nx = [nx x]; ny = [ ny y ];
  fprintf(1,'[x,y] = [%g %g]\n',x,y);
end
