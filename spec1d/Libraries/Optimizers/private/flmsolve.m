function [xf, S, nfJ,exitflag, J] = flmsolve(varargin)
%FLMSOLVE Summary of this function goes here
%   Detailed explanation goes here

%   LMFsolve(FUN,Xo,fixed,Options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUN = varargin{1};            %   function handle
if ~(isvarname(FUN) || isa(FUN,'function_handle'))
    error('FUN Must be a Function Handle or M-file Name.')
end

xc    = varargin{2};             %   Xo
fixed = varargin{3};

if nargin > 3                 %   OPTIONS
    if isstruct(varargin{4})
        options = varargin{4};
    else
        if ~exist('options','var')
            options = flmsolve_fast('default');
        end
        for i = 4:2:size(varargin,2)-1
            options = flmsolve_fast(options, varargin{i},varargin{i+1});
        end
    end
else
    if ~exist('options','var')
        options = flmsolve_fast('default');
    end
end

x  = xc(:);
lx = length(x);

r = sqrt(feval(FUN,x));             % Residuals at starting point
if length(r) == 1, 
    r = r*ones(5,1)/5;
else
    r=r(:); 
end
%~~~~~~~~~~~~~~

S = r'*r;
epsx = options.XTol(:);
J = options.Jacobian(FUN,r,x,epsx);
nfJ = lx+1;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A = J.'*J;                    % System matrix
v = J.'*r;

D = options.ScaleD;
if isempty(D)
    D = diag(diag(A));        % automatic scaling
    for i = 1:lx
        if D(i,i) == 0
            D(i,i)=1; 
        end
    end
else
    if numel(D) > 1
        D = diag(sqrt(abs(D(1:lx)))); % vector of individual scaling
    else
        D = sqrt(abs(D))*eye(lx);     % scalar of unique scaling
    end
end

Rlo = 0.25; 
Rhi = 0.75;
l   = 1;
lc  = 0.75;     
cnt = 0;
% d = options.XTol;             %   vector for the first cycle

exitflag = 0;

while exitflag %   MAIN ITERATION CYCLE
    try
        d = pinv(A+l*D)*v;
    catch
        d = (A+l*D)\v;            % negative solution increment
    end
    d(~fixed) = zeros(1,sum(~fixed));
    xd = x-d;
    rd = sqrt(feval(FUN,xd));
    if length(rd) == 1, 
        rd = rd*ones(5,1)/5;
    else
        rd = rd(:); 
    end
    nfJ = nfJ + 1;
    %   ~~~~~~~~~~~~~~~~~
    Sd = rd.'*rd;
    dS = d.'*(2*v - A*d);       % predicted reduction
    R = (S - Sd)/dS;
    
    if R > Rhi
        l = l/2;
        if l < lc
            l = 0; 
        end
    elseif R < Rlo
        nu = (Sd - S)/(d.'*v) + 2;
        if nu < 2
            nu = 2;
        elseif nu > 10
            nu = 10;
        end
        if l == 0
            lc = 1/max(abs(diag(pinv(A))));
            l = lc;
            nu = nu/2;
        end
        l = nu*l;
    end
    cnt = cnt+1;
    
    S = Sd; 
    x = xd; 
    r = rd;
    J = options.Jacobian(FUN,r,x,epsx);
    J(:,~fixed) = 0;
    nfJ = nfJ+lx;
    %       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A = J.'*J;   
    v = J.'*r;
    
    
    if      cnt >= options.MaxIter
        exitflag = -2; % max iteration reached
    elseif  all(abs(d) < options.XTol/1000)   
        exitflag = -1; % parameter change increment is negligible
    elseif  all(abs(r) < options.FunTol) 
        exitflag = -5; % function change increment reached
    elseif  nfJ > options.Evals
        exitflag = -3; % max nb function evaluations reached
    end
end
xf = x;                         %   final solution

end

