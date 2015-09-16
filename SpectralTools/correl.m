function [y, D, R] =  correl( d, r, n)
% correl : Correlation
% Syntax: [y, fft1, fft2] = correl ( data1, data2, {n=[]} ) or correl(data)
%
% Computes correlation function of 'data1' and 'data2', or autocorrelation.
% Side effects and periodicity artifacts are managed by zero padding.
% If n is supplied by user, a n-points FFT is used, else using powers of 2.
% FFT's are also returned.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: correlation

% uses : fft, ifft
% From : Numerical Recipes in C. p545
% Part of 'Spectral tools'. E.Farhi. 07/96

if (nargin < 1)
	usage('y = correl ( data1, data2, {n=[]} ) or correl(data)');
end

if (nargin == 1)
	r=d;
end

if min(size(d) ~= 1) | min(size(r) ~= 1)
	error('Should be used with vectors.');
end

lr = length(r);
ld = length(d);

if ~exist('n') 
	n = [];
end

if isempty(n)	% just higher power of 2
	n = 2^(ceil(log (lr + ld -1) / log(2)));
end

D = fft(d, n);
R = fft(r, n);

% -------------------- computes -----------------


y = D.*conj(R);		% correlation

% -------------------- end computes -----------------

y = ifft(y);

y = y(1:ld);

% Final cleanups:  if both x and b are real respectively integer, y
% should also be

  if (~ (any (imag (d)) | any (imag (r))))
    y = real (y);
  end
  if (~ (any (d - round (d)) | any (r - round (r))))
    y = round (y);
  end


	
