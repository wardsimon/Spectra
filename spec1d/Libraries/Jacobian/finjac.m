%   FINJAC       numerical approximation to Jacobi matrix
%   %%%%%%
function J = finjac(FUN,r,x,epsx)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pars=column, function=row vector or scalar
lx = length(x);
J  = zeros(length(r), lx);
x   = x(:)';
r   = r(:);

if numel(epsx) == 1, 
    epsx = epsx*max(abs(x),1); 
end
if any(epsx == 0)
    epsx(~epsx) = 1e-4; 
end

epsx = epsx(:)';

for k = 1:lx
    dx = .5*epsx(k); % !!FIX!! was 0.25*epsx(k)
    xd = x;
    xd(k) = xd(k) + dx;
    rd    = sqrt(feval(FUN, xd));
    if length(rd) == 1
        rd = rd*ones(5,1)/5;
    else
        rd = rd(:);
    end
    %   ~~~~~~~~~~~~~~~~
    J(:,k) = ((rd-r)/dx);
end
end % finjac