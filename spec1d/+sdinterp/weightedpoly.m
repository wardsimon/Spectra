function [p, S] = weightedpoly(x,y,e,n)
% WEIGHTEDPOLY 1-D polynomial interpolation weighted by error
% Simon Ward 01/02/2016
c = onCleanup(@() warning('on','MATLAB:nearlySingularMatrix'));

x = x(:);
y = y(:);
w = 1./e(:);

% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
    V(:,j) = x.*V(:,j+1);
end
% Solve least squares problem.
if verLessThan('matlab','8.3')
    [Q, R] = qr_lite(bsxfun(@times,V,w));
else
    [Q,R] = qr(bsxfun(@times,V,w),0);
end

warning('off','MATLAB:nearlySingularMatrix')
p = R\(Q'*bsxfun(@times,y,w));    % Same as p = V\y;
warning('on','MATLAB:nearlySingularMatrix')

if nargout > 1
    r = y - V*p;
    % S is a structure containing three elements: the triangular factor from a
    % QR decomposition of the Vandermonde matrix, the degrees of freedom and
    % the norm of the residuals.
    S.R = R;
    S.df = max(0,length(y) - (n+1));
    S.normr = norm(r);
end

end

function [Q, R] = qr_lite(A)
%        Reference:
%        N. J. Higham, Accuracy and Stability of Numerical Algorithms,
%        Second edition, Society for Industrial and Applied Mathematics,
%        Philadelphia, PA, 2002; sec 19.8.
[m, n] = size(A);
Q = zeros(m,n);
R = zeros(n);

for k=1:n
    R(k,k) = norm(A(:,k));
    Q(:,k) = A(:,k)/R(k,k);
    R(k,k+1:n) = Q(:,k)'*A(:,k+1:n);
    A(:,k+1:n) = A(:,k+1:n) - Q(:,k)*R(k,k+1:n);
end
end