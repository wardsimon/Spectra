function v=vect2row(vin)
% vect2row : Make row
%Syntax: V=vect2row(Vin)
%
% 'V' is 'Vin' reshaped as a row vector.
% For matrix input, output is as Nrows < Ncolumns.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: make row

% E.Farhi.   12/95   (manuf@ldv.univ-montp2.fr)
% Spectral tools.

[nr,nc] = size (vin);

if ( nr > nc)
	v=vin';
else
	v=vin;
end

%end
