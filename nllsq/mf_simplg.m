function [p, sig,f,v,fv,covp,corp,r2,rv]=mf_simplx(x,y,err,pin,notfixed,func,fcp)
%
% MFIT simplex routine (see fmins matlab function)
% FMINS uses a Nelder-Mead type simplex search method.
% syntax : [p, sig, f]=mf_simplx(x,y,err,pin,notfixed,func,{fcp})
% with fcp = [ dp niter stol ]
% When fcp is given, outputs are disabled (silent mode).

% This version adds an option to the data window to
% allow the fit control parameters to be changed.
% EF 03.98

flagout = 1;
showgraph = 1;
if nargin < 7, fcp = []; end

%----- Add fit control to data window if necessary

[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;
hmf_fcpmenu=findobj('Tag','mf_fcpmenu');

if hmf_data~=0 & isempty(fcp)
  if isempty(hmf_fcpmenu) % create menu
     figure(hmf_data);
     hmf_fcpmenu=uimenu(hmf_data,'Label','Control','Tag','mf_fcpmenu');
     uimenu(hmf_fcpmenu,'Label','Change',...
            'Callback','mf_fcp');
     fcp=[0.005 20 0.00001];
     set(hmf_fcpmenu,'Userdata',fcp);
  else
     fcp=get(hmf_fcpmenu,'Userdata');
  end
else
  flagout = 0;
  showgraph = 0;
end

if isempty(fcp)
  fcp=[0.1 20 0.1];
end

% define starting values ------------------------------

if length(notfixed) ~= length(pin)
  notfixed = ones(length(pin),1);
end

if length(x) ~= length(y)
  x=1:length(y);
end

if length(err) ~= length(y)
  err = ones(length(x),1);
end

selected = find(err);

wt=0*y;
wt(selected)=1./err(selected);
x=x(:);y=y(:); wt=wt(:); pin=pin(:); notfixed = notfixed(:);
m=length(y); n=length(pin); p = pin;
dp=notfixed*fcp(1); niter=fcp(2); stol=fcp(3);
i = find(notfixed & (abs(dp) <= 1e-12));  % for small not fixed params
if ~isempty(i)
  dp(i) = max(2e-12,abs(min(pin(find(pin)))*fcp(1)/1000));
end

fixed = find(~notfixed);

sig = []; covp = []; corp = []; r2 = []; rv = []; f = []; v = []; fv = [];
if isempty(pin)
  return
end

if (m <= n)
  if flagout
    disp('Warning : need to extend data...');
  end
  newm = ceil(n/m+1)*m;
  newx = linspace(x(1), x(length(x)), newm);
  newx = newx(:);
  y = interp1(x,y,newx);
  wt = interp1(x,wt,newx);
  wt = wt(:);
  x=newx;
  m = newm;
end

% interface/call to fmins --------------------------------------------------

funfcn=func;
xsav=x;
ysav=y;

x=pin;
fixed=find(dp == 0);
notfixed = (dp ~= 0);
tol=dp;
tol2=stol * length(ysav);  %  * length(ysav);

% f=feval(funfcn,xsav,x); r=wt.*(ysav-f); f=r'*r;

if flagout
  disp('* Nelder-Mead type simplex fit')
end

% FMINS  --------------------------------------------------

n = prod(size(x));

% Set up a simplex near the initial guess.
xin = x(:); % Force xin to be a column vector
v = xin;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
x(:) = v;
f=feval(funfcn,xsav,x); f=f(:); r=wt.*(ysav-f); fv=r'*r;
pbest=x;
ssbest=fv;
fbest=f;
niter=niter*length(find(dp));
if flagout & showgraph
  figure(hmf_data);
  hfit=line(xsav,f,'erasemode','xor','color','c','Tag','mf_fitline');
end

tic

% Following improvement suggested by L.Pfeffer at Stanford
usual_delta = tol;             % percent deltas for non-zero terms
zero_term_delta = min(tol(find(tol)))/50;      % Even smaller delta for zero elements of x
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

   y = xin;
   if y(j) ~= 0
        y(j) = (1 + 2*(rand-0.5)*usual_delta(j))*y(j);
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
if flagout > 0
  disp(sprintf('* Beginning fit (max %d iterations)',niter));
  disp('--------------------------------------')
  disp('Iteration  Time(s)  Residual Variance');
  fprintf(1,'   %3d      %6.2f   %8.3f %s\n', func_evals, toc, fv(1)/length(ysav),'initial  (best)');
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
    f = feval(funfcn,xsav,x); f=f(:); r=wt.*(ysav-f); fr=r'*r;
    func_evals = func_evals + 1;
    vk = vr;  fk = fr; how = 'reflect ';
    if fr < fv(n)
        if fr < fv(1)
            ve = gamma*vr + (1-gamma)*vbar;
            ve(fixed)=pin(fixed);
            x(:) = ve;
            f=feval(funfcn,xsav,x); f=f(:); r=wt.*(ysav-f); fe=r'*r;
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
        f=feval(funfcn,xsav,x); f=f(:); r=wt.*(ysav-f); fc=r'*r;
        func_evals = func_evals + 1;
        if fc < fv(n)
            vk = vc; fk = fc;
            how = 'contract';
        else
            for j = 2:n
                v(:,j) = (v(:,1) + v(:,j))/2;
                x(:) = v(:,j);
                f=feval(funfcn,xsav,x); f=f(:); r=wt.*(ysav-f); fv(j)=r'*r;
            end
        func_evals = func_evals + n-1;
        vk = (v(:,1) + v(:,n+1))/2;
        vk(fixed)=pin(fixed);
        x(:) = vk;
        f=feval(funfcn,xsav,x); f=f(:); r=wt.*(ysav-f); fk=r'*r;
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
    if (flagout) fprintf(1,'   %3d      %6.2f   %8.3g %s', iter, toc, ss,how); end
    if (ss <= ssbest) pbest=x; ssbest=ss; fbest=f;
  if (flagout)
    if showgraph
      mf_upars(x,sig);
      set(hfit,'Ydata',f);
      drawnow;
    end
    fprintf(1,' (best)');
  end
    end
    if (flagout) fprintf(1,'\n'); end

end
x(:) = v(:,1);
if func_evals>=niter

end
if func_evals==niter
  if flagout
    disp(' Getting best fit... ')
  end
  kvg=0;
else
  kvg=1;
end

% --------------------------------------------------

% set return values
%
p=pbest(:);
f=fbest;

y=ysav;
x=xsav;

r=corrcoef(y.*wt,f.*wt);
r2=r.*r';
ss=sum(((f-ysav).*wt).^2/length(ysav));

% now evaluates 'v' params dispersion near final solution : statistic evaluation in simplex
% v is nparams columns
if nargout > 4
  [sig,covp,corp,r2,rv,f] = mf_sig(x,y,wt,p,p*0.01,notfixed,func,f);
else
  sig2=std(v');
  sig = sig2(:).*notfixed(:);
end
fv = fv/length(y); % get residual variances
if flagout
  fprintf(1,'Best RV : %.3f - Worse RV : %.3f \n', min(fv),max(fv));
end


%-------- Work out Chi squared ------------------------
v=length(y)-length(find(dp));   % # of degrees of freedom (#sel points - #free pars)
ChiSq = sum( ((y-f).*wt).^2 )/v;




