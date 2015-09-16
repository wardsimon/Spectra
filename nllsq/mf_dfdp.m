function y = dfdp(x,f,p,dp,func)
%df/dp = dfdp(x,y,p,dp,func)
%
% Returns the partial derivatives of function 'func'.
% 'x'(vect) is x axis values, 'y' is y values, 'p' and 'dp' are parameters and
% their variation. 'func' is the function concerned (y=func(x,p)).
% output 'df/dp' is a vector(or matrix) of partials varying of 'dp'.

% uses : none
% Part of 'Spectral tools'.
% *** Do not modify *** or use a copy. ***
% E.Farhi.   12/95  
% You can either enter exact expression of df/dp, or compute it
% with finite differences.

% example of direct computation for func : y=p(1)*exp(-p(2)*x)
%y=[exp(-p(2)*x) -p(1)*x.*exp(-p(2)*x)];
%return;

% finite differencies.

x2 = x;

[nr,nc] = size (x);
if (nr < nc)
	x=x';
end

p0 = p; 		% save init params
y = zeros(length(x), length(p));
i = find((abs(dp) <= 1e-12) | isnan(dp) | isinf(dp));	% thoses params should not vary
if ~isempty(i)
	dp(i) = 0;
end

for i=1:length(p)
	if (dp(i) ~= 0)
		p(i) = p(i) + dp(i);
		t=feval(func,x,p);
		if length(t) ~= length(f)
			disp('Warn : fix length problem in jacobian')
			t = t(1:length(x));
		end
		y(:,i) = (t(:)-f(:))/dp(i);
		p=p0;
	end
end

if (nr < nc)
	y = y';
end

