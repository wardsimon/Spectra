function a=is_row(x)
% is_row : true when row vector
%Syntax: a=is_column(x)
%
%True if 'x' is a column vector.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: true when row vector

[nx,ny]=size(x);
a=(ny >= nx);
%end
