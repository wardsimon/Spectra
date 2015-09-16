function [ny,c] = mf_smooth(yd,N,M)
%Data smoothing by Savitzky-Golay method
%Syntax: [smoothed_y, coefs] = smooth(y,{N,M})
%
% Smoothes the y signal by M-th order Savitzky-Golay method with N points.
% This algorithm assumes that corresponding x axis is evenly spaced.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  data smoothing by Savitzky-Golay method

% Part of 'Spectral tools'. E.Farhi. 07/96
% From : Numerical recipes in C. p 650

% Argument processing -----------------------------------------------

if (nargin < 1)
	error('usage: [smoothed_y, coefs] = smooth(y,{N,M=2})');
end

nd=length(yd);

if (nargin < 3)
	M = 2;
end

if (nargin <= 1)
	N = ceil(max(2*M,nd/50));
end

if (nd <= 2*N)
	error('not enough points in data');
end

%  Savitzky-Golay coefficients -----------------------------------------------

N = ceil(N/2);

A = zeros (2*N +1, M+1);
c = zeros(1,2*N+1);

n=(-N):N;
for j=0:M
    A(:,j+1) = n'.^j;		% Aij = i^j
end

B = pinv(A'*A);
B = B(1,1:(M+1));

for n=1:(2*N+1)
  c(n) = A(n,:) * B';	% these are Savitzky-Golay coefficients
end

% Smoothing ------------------------------------------------------------------

ny = 0*yd;

for n=1:(2*N+1)
  ny((N+1):(nd-N)) = ny((N+1):(nd-N)) + c(n) * yd(n:(nd-2*N-1+n));
end
for n=1:N
  ny(n)=yd(n);
  ny(nd-n+1) = yd(nd-n+1);
end

ny = ny(:);
