function [p, sig, f, covp,corp,r2,rv]=mf_flsqr(x,y,err,pin,notfixed,func,fcp)
% mf_flsqr : Marquart Levenberg fit routine (fast, non graphic)
% Version 3.beta
% Levenberg-Marquardt nonlinear regression of f(x,p) to y(x)
% syntax : [p, sig, f]=mf_flsqr(x,y,err,pin,notfixed,func,{fcp})
% with fcp = [ dp niter stol ]
% When fcp is given, outputs are disabled (silent mode).

% Richard I. Shrager (301)-496-1122
% Modified by A.Jutan (519)-679-2111
% Modified by Ray Muzic 14-Jul-1992
% DFM 24.3.96 % EF 03.98

% 0 error means point is not selected.

flagout = 1;
showgraph = 0;
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
  fcp=[0.01 20 0.001];
end

% define starting values ------------------------------

if length(notfixed) ~= length(pin)
  notfixed = ones(length(pin),1);
end

if length(x) ~= length(y)
  x=1:length(y);
end

if length(err) ~= length(y)
  err = ones(length(y),1);
end

selected = find(err);

wt=0*y;
wt(selected)=1./err(selected);
x=x(:);y=y(:); wt=wt(:); pin=pin(:); notfixed = notfixed(:);
m=length(y); n=length(pin);

dp=abs(notfixed*fcp(1).*pin); niter=fcp(2); stol=fcp(3);
i = find(notfixed & (abs(dp) <= 1e-12));  % for small not fixed params
if ~isempty(i)
  dp(i) = max(2e-12,abs(min(pin(find(pin)))*fcp(1)/1000));
end

sig = []; covp = []; corp = []; r2 = []; rv = []; f = [];
if isempty(pin)
  return
end

% if m <= n : SVD method can't work -> we extend data by interpolation

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

ldfdp = 'mf_dfdp';
p = pin;

% get into leasqr.m standard routine ------------------------------

options=[zeros(n,1) Inf*ones(n,1)];
nor = n; noc = 2;

pprec=options(:,1);
maxstep=options(:,2);

% set up for iterations (with initial data) -------------------------

f=feval(func,x,p); f=f(:); wt = wt/sum(wt);
fbest=f; pbest=p;   % initial function values
r=wt.*(y-f);          % error
sbest=r'*r;         % square error
ss = sbest;
nrm=zeros(n,1);
chgprev=Inf*ones(n,1);
kvg=0;
epsLlast=1;
epstab=[.1 1 1e2 1e4 1e6 ];

% do iterations -----------------------------------------------------

tic

if flagout >0
  disp('* Levenberg-Marquardt nonlinear regression')
  if showgraph
    [hmf_ctrl, hmf_data, hmf_pars]=mf_figs;
    figure(hmf_data);
    hfit=line(x,f,'erasemode','xor','color','c','Tag','mf_fitline');
  end
  disp(sprintf('* Beginning fit (max %d iterations)',niter));
  disp('--------------------------------------')
  disp('Iteration  Time(s)  Residual Variance');
end

tic
for iter=1:niter,
  pprev=pbest;          % current params
  prt=feval(ldfdp,x,fbest,pprev,dp,func); % partials of func for x.
%  prt = prt(:);
  r=wt.*(y-fbest);        % error
  sprev=sbest;          % square error
  sgoal=(1-stol)*sprev;       % sq error to achieve
  for j=1:n,      % computes partials norm
    if dp(j)==0,
      nrm(j)=0;
    else
      prt(:,j)=wt.*prt(:,j);
      nrm(j)=prt(:,j)'*prt(:,j);
      if nrm(j)>0,
        nrm(j)=1/sqrt(nrm(j));
      end
    end
    prt(:,j)=nrm(j)*prt(:,j); % normalizes partials
  end
  [prt,s,v]=svd(prt,0);   % prt=unit matr, s=eigenval : prt=prt*s*v'.
  s=diag(s);      % get eigenvalues (diagonal of sigma)
  g=prt'*r;     % gradient by Gauss-Newton formula.
  for jjj=1:length(epstab),
    epsL = max(epsLlast*epstab(jjj),1e-7);
    se=sqrt((s.*s)+epsL);
    gse=g./se;
    chg=((v*gse).*nrm);   % change on params
        % check the change constraints and apply as necessary
    ochg=chg;
    for iii=1:n,
      if (maxstep(iii)==Inf), break; end;
      chg(iii)=max(chg(iii),-abs(maxstep(iii)*pprev(iii)));
      chg(iii)=min(chg(iii),abs(maxstep(iii)*pprev(iii)));
    end;
    if (any(ochg ~= chg)),
      if (flagout >0)
        disp(['Change in parameter(s): ' ...
         sprintf(1,'%d ',find(ochg ~= chg)) 'were constrained']);
      end
    end;
    aprec=abs(pprec.*pbest);       %---
    if (any(abs(chg) > 0.1*aprec)) %---  % only worth evaluating function if
      p=chg+pprev;                        % there is some non-miniscule change
      f=feval(func,x,p); f=f(:);
      if showgraph
  set(hfit,'Ydata',f);
  drawnow;
      end
      r=wt.*(y-f);
      ss=r'*r;
      if ss<sbest,
        pbest=p;
        fbest=f;
        sbest=ss;
      end;
      if ss<=sgoal,
        break;
      end;
    end;                          %---
  end;
  epsLlast = epsL;
  if ss<eps,
    break;
  end
  aprec=abs(pprec.*pbest);
%  [aprec chg chgprev]
  if (all(abs(chg) < aprec) & all(abs(chgprev) < aprec)),
    kvg=1;
    if (flagout)
      fprintf(1,'Parameter changes converged to specified precision\n');
    end
    break;
  else
    chgprev=chg;
  end;
  if ss>sgoal,
    break;
  end;
  if flagout
  mf_upars(p,[]);
    disp(sprintf('   %3d      %6.2f   %8.3g', iter, toc, ss/length(y)));
  end
end

% set return values
%
p=pbest;
f=fbest;
ss=sbest;
kvg=((sbest>sgoal)|(sbest<=eps)|kvg);
if kvg ~= 1 & flagout, disp(' CONVERGENCE NOT ACHIEVED! '), end;

if nargout > 1
  [sig,covp,corp,r2,rv,f] = mf_sig(x,y,wt,p,p*0.01,notfixed,func,f);
  disp('* Covariance matrix of estimated parameters')
  disp(covp)
  disp('* Correlation matrix of estimated parameters ')
  disp(corp)
end
