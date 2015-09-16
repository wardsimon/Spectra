function [pout,deltapout,r2,ss,f,func_evals,kvg,ChiSq,v,fv]=simplex(x,ys,pin,func,wt,dp,niter,stol,last)
% simplex : Simplex routine (see fmins matlab function)
% FMINS uses a Nelder-Mead type simplex search method.
% Syntax: [p,deltap,r2,rv,f,iter,kvg,ChiSq,v,fv]=
%                   simplex(x,y,Pin,Func,{wt,dp,niter,stol,mode='normal'})
%
% INPUT:
%  x,y (vect)    : data to be fitted.
%  Pin (vect)    : initial values for parameter search.
%  Func (string) : function name (between quotes), of the form 'f(x,p)'.
%
% OUTPUT:
%  p (vect)      : final parameters.
%  r2 (scalar)   : correlation coefficient for fit.
%  deltap (v)    : sigma on parameters
%  rv (scalar)   : residual variance.
%  f (vect)      : function values at the end of fitting process.
%  ChiSq         : Chi Squarred
%
%                            ---------------
% OPTIONS IN:
%  wt (s or v)   : statistical weights for error calculation (1/sqrt(y) is defualt).
%  dp (s or v)   : variation window for parameters, 0 for fixed param (0.1 is default, 0 to fix param)
%  niter (integ) : maximal number of iterations (20 is deualt).
%  stol  (scal)  : tolerency on square error sum fractional improvement (0.1 is default).
%  mode          : 'silent', 'min', 'normal' or 'verbose' (do plot)
%
% OPTIONS OUT:
%  iter (integer): number of iterations processed.
%  kvg (boolean) : convergence flag (1 if ok).
%  v (matrix)    : set of parameters used in fit
%  fv (vector)   : set of RV for each set in 'v'

% EF 09.97 from numerical recipies in C and FMINS matlab function

% argument processing -----------------------------------------------

if nargin < 9, last = []; end
if nargin < 8, stol = []; end
if nargin < 7, niter = []; end
if nargin < 6, dp = []; end
if nargin < 5, wt = []; end


if isempty(last)
  last = 'normal';
end

if strcmp(last,'silent')
  tmp = 0;
elseif strcmp(last,'min')
  tmp = 1;
elseif strcmp(last,'verbose')
  tmp = 3;
else
  tmp = 2;
end

if tmp>0
  fprintf(1,'\nSimplex Multivariable Fit\n');
end


if isempty(dp)
  dp=.1;
  if (tmp>1)
                            disp('dp = 0.1 (10 %) ');
  end
end
if isempty(wt)
  wt=1./sqrt(abs(ys+(ys == 0)));
  if (tmp>1)
    disp('using standard normalisation (1/sqrt(y))')
  end
end           % data weights
if isempty(niter)
  niter=20;
  if (tmp>1)
    fprintf(1,'Iterations : %i\n',niter);
  end
end           % max nb of iters
if isempty(stol)
  stol=.1;
  if (tmp>1)
    fprintf(1,' Tolerancy : %g\n',stol);
  end
end           % tolerancy

if length(wt) == 1
  wt = wt * ones(length(ys),1);
end

nzerow = find(wt);
zerow = find(wt == 0);
if ~isempty(zerow)
  if tmp>1
    disp('Get selected points');
  end
  x = x(nzerow);
  ys = ys(nzerow);
  wt = wt(nzerow);
end

%change all vectors to columns
x=x(:);
y=ys(:);
wt=wt(:);
pin=pin(:);
dp=dp(:);
% check data vectors- same length?
m=length(y); n=length(pin); p=pin;[m1,m2]=size(x);
if (m1~=m)
  error('input(x)/output(y) data must have same number of rows ');
end

if (m > 10000)
    warning('HUH ! That makes a lot of points ! I may crash during fit...');
end

if (m <= n)
        disp('Warning : need to extend data...');
        newm = ceil(n/m+1)*m;
        newx = linspace(x(1), x(length(x)), newm);
        newx = newx(:);
        y = interpsp(x,y,newx);
        wt = interpsp(x,wt,newx);
        x=newx;
        m = newm;
end

if (length(dp) == 1)
  dp=dp*ones(n,1);
end

% interface/call to fmins --------------------------------------------------

funfcn=func;
xsav=x;
ysav=y;

x=pin;
fixed=find(dp == 0);
notfixed = (dp ~= 0);
tol=dp;
tol2=stol;

% f=feval(funfcn,xsav,x); r=wt.*(ysav-f); f=r'*r;

% FMINS  --------------------------------------------------

n = prod(size(x));

% Set up a simplex near the initial guess.
xin = vect2column(x); % Force xin to be a column vector
v = xin;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
x(:) = v;
f=vect2column(feval(funfcn,xsav,x)); r=wt.*(ysav-f); fv=r'*r;
pbest=x;
ssbest=fv;
fbest=f;
niter=niter*length(find(dp));

tic

% Following improvement suggested by L.Pfeffer at Stanford
usual_delta = tol;             % percent deltas for non-zero terms
zero_term_delta = tol/50;      % Even smaller delta for zero elements of x
for j = find(dp')
   y = xin;
   if y(j) ~= 0
        y(j) = (1 + usual_delta(j))*y(j);
   else
        y(j) = zero_term_delta;
   end
%   y(fixed)=pin(fixed);
   v = [v y];
   x(:) = y; f = feval(funfcn,xsav,x); r=wt.*(ysav-f); f=r'*r;
   fv = [fv  f];

   y = xin;
   if y(j) ~= 0
        y(j) = (1 - usual_delta(j))*y(j);
   else
        y(j) = -zero_term_delta;
   end
%   y(fixed)=pin(fixed);
   v = [v y];
   x(:) = y; f = feval(funfcn,xsav,x); r=wt.*(ysav-f); f=r'*r;
   fv = [fv  f];
end
[fv,j] = sort(fv);
v = v(:,j);

[np1,n] = size(v);
n = n-1;
func_evals = n+1;
if tmp > 0
disp('')
disp(sprintf('*Beginning fit (max %d iterations)',niter));
if tmp > 1
disp('--------------------------------------')
disp('Iteration  Time(s)  Residual Variance');
fprintf(1,'   %3d      %6.2f   %8.3f %s\n', func_evals, toc, fv(1)/length(ysav),'initial  (best)');
end
end

alpha = 1;  beta = 1/2;  gamma = 2;
onesn = ones(1,n);
ot = 2:n+1;
on = 1:n;

% Iterate until the diameter of the simplex is less than tol.
while func_evals < niter
    if max(abs(fv(1)-fv(ot))) <= tol2
        break
    end

    % One step of the Nelder-Mead simplex algorithm

    vbar = (sum(v(:,on)')/n)';
    vr = (1 + alpha)*vbar - alpha*v(:,n+1);
    vr(fixed)=pin(fixed);
    x(:) = vr;
    f = vect2column(feval(funfcn,xsav,x)); r=wt.*(ysav-f); fr=r'*r;
    func_evals = func_evals + 1;
    vk = vr;  fk = fr; how = 'reflect ';
    if fr < fv(n)
        if fr < fv(1)
            ve = gamma*vr + (1-gamma)*vbar;
            ve(fixed)=pin(fixed);
            x(:) = ve;
            f=vect2column(feval(funfcn,xsav,x)); r=wt.*(ysav-f); fe=r'*r;
            func_evals = func_evals + 1;
            if fe < fv(1)
                vk = ve; fk = fe;
                how = 'expand  ';
            end
        end
    else
        vt = v(:,n+1); ft = fv(n+1);
        if fr < ft
            vt = vr; ft = fr;
        end
        vc = beta*vt + (1-beta)*vbar;
        vc(fixed)=pin(fixed);
        x(:) = vc;
        f=vect2column(feval(funfcn,xsav,x)); r=wt.*(ysav-f); fc=r'*r;
        func_evals = func_evals + 1;
        if fc < fv(n)
            vk = vc; fk = fc;
            how = 'contract';
        else
            for j = 2:n
                v(:,j) = (v(:,1) + v(:,j))/2;
                x(:) = v(:,j);
                f=vect2column(feval(funfcn,xsav,x)); r=wt.*(ysav-f); fv(j)=r'*r;
            end
        func_evals = func_evals + n-1;
        vk = (v(:,1) + v(:,n+1))/2;
        vk(fixed)=pin(fixed);
        x(:) = vk;
        f=vect2column(feval(funfcn,xsav,x)); r=wt.*(ysav-f); fk=r'*r;
        func_evals = func_evals + 1;
        how = 'shrink  ';
        end
    end
    v(:,n+1) = vk;
    fv(n+1) = fk;
    [fv,j] = sort(fv);
    v = v(:,j);

    iter=func_evals;
    ss=sum(((f-ysav).*wt).^2/length(ysav));
    if (tmp > 1) fprintf(1,'   %3d      %6.2f   %8.3f %s', iter, toc, ss,how); end
    if (ss <= ssbest) pbest=x; ssbest=ss; fbest=f;
  if (tmp>1) fprintf(1,' (best)'); end
    end
    if (tmp>1) fprintf(1,'\n'); end

end
x(:) = v(:,1);
if func_evals>=niter

end
if func_evals==niter
  if tmp>0
    disp(' Getting best fit... ')
  end
  kvg=0;
else
  kvg=1;
end

% --------------------------------------------------

% set return values
%
p=pbest;
f=fbest;

y=ysav;
x=xsav;

r=corrcoef(y.*wt,f.*wt);
r2=r.*r';
ss=sum(((f-ysav).*wt).^2/length(ysav));

% now evaluates 'v' params dispersion near final solution : statistic evaluation in simplex
% v is nparams columns

sig2=std(v');
sig = sig2.*notfixed';
fv = fv/length(y); % get residual variances
if tmp>0
  fprintf(1,'Best RV : %.3f - Worse RV : %.3f \n', min(fv),max(fv));
end


%-------- Work out Chi squared ------------------------
v=length(y)-length(find(dp));   % # of degrees of freedom (#sel points - #free pars)
ChiSq = sum( ((y-f).*wt).^2 )/v;

pout = vect2column(p);
deltapout = vect2column(sig);

if (tmp>0)
  disp( '* Least Square Estimates of Parameters')
  disp('      value      sigma(+/-)')
  disp([ pout deltapout ])
  disp( '* Correlation Coefficient for fit R^2')
  disp(r2)
  disp( '* Residual Variance')
  disp(ss);
  disp('* Chi squared')
  disp(ChiSq)
end

if (tmp>2)
  fprintf(1,'Fit results sent to graph...\n');
  title('Fit result');
  xlabel('X'); ylabel('Y');
  plot(x(:,1),y,x(:,1),f);  % plot command
end

p=p';
if (size(ys,1)==1)
  f = f';
end



