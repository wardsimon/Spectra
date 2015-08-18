function a=is_column(x)
% is_column : true when column vector
%Syntax: a=is_column(x)
%
%True if 'x' is a column vector.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: true when column vector

[nx,ny]=size(x);
a=(ny <= nx);
%end
