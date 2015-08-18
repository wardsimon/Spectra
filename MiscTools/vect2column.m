function v=vect2column(vin)
% vect2column : Make column
%Syntax: V=vect2column(Vin)
%
% 'V' is 'Vin' reshaped as a column vector.
% For matrix input, output is as Nrows > Ncolumns.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: make column

% E.Farhi.   12/95   (manuf@lsmv.univ-montp2.fr)
% Spectral tools.

[nr,nc] = size (vin);

if ( nr < nc)
	v=vin';
else
	v=vin;
end;

%end
