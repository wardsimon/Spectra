function [new_x, p] = axisrescale(x, pmat, order,last)
% axisrescale : Rescales an axis by Lagrange interpolation
%Syntax: [new_x, p] = axisrescale(x, pmat, {order=2,mode='silent'})
%
% This function rescales an x axis using polynomial interpolation.
% 'pmat' is a matrix of rows [old_x new_x].
% returns new axis and interpolation polynome.

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description: rescales an axis by Lagrange interpolation

% uses : interp.m
% Part of 'Spectral tools'. E.Farhi. 07/96

if (nargin < 2)
	error('usage : new_x = axisrescale(x, pmat, {order=2,mode=''silent''})');
end

if ~exist('order')
	order = 2;
end

if isempty(order)
	order = 2;
end

if ~exist('last')
	last = 'silent';
end

if  isempty(last)
	last = 'silent';
end

if strcmp(last,'silent')
	tmp = 0;
else
	tmp = 1;
end

np = size(pmat,1);
if (np < 2)
	error('need at least 2 points in interpolation matrix.');
end
lx = length(x);
if tmp
	fprintf(1,'\nusing order %i rescaling on %i points.\n', order, lx);
end

[a, p] = interpsp(pmat(:,1) , pmat(:,2), pmat(ceil(np/2),1), np, order);

new_x = polyval(p, x);
if tmp
	fprintf(1,'interpolation polynome is [ n..1 ] :');
	disp(p);
end
