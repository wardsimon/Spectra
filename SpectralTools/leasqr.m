function [pout,deltapout,r2,covp,corp,ss,f,iter,kvg,ChiSq]=leasqr(x,ys,pin,func,wt,dp,niter,stol,ldfdp,options,last)
% leasqr : Non Linear Least Square multivariable fit.
%Syntax: [p,deltap,r2,covp,corp,rv,f,iter,kvg,ChiSq]=
%                   leasqr(x,y,Pin,Func,{wt,dp,niter,stol,dfdp,options,mode='normal'})
%
% Levenberg-Marquardt nonlinear regression of func(x,p) to y(x).
%
% INPUT:
%  x,y (vect)    : data to be fitted.
%  Pin (vect)    : initial values for parameter search.
%  Func (string) : function name (between quotes), of the form 'f(x,p)'.
%
% OUTPUT:
%  p (vect)      : final parameters.
%  deltap (v)    : sigma on parameters = sqrt(diag(covp))
%  r2 (scalar)   : correlation coefficient for fit.
%  covp (matrix) : error and crossed error on parameters.
%  corp (matrix) : correlation between parameters.
%  rv (scalar)   : residual variance.
%  f (vect)      : function values at the end of fitting process.
%  ChiSq         : Chi Squarred
%  
%                            ---------------
% OPTIONS IN:
%  wt (s or v)   : statistical weights for error calculation (1/sqrt(y) is default).
%  dp (s or v)   : fractional incr of params for numerical partials (0.01 is default, 0 to fix param).
%  niter (integ) : maximal number of iterations (20 is default).
%  stol  (scal)  : tolerency on square error sum fractional improvement (1e-3 is default).
%  dfdp (string) : partials function name of Func ('dfdp' is default).
%  options (matr): col1 = desired fractional precision in parameter estimates.
%                  col2 = maximum fractional step change in parameter vector.
%  mode          : 'silent', 'min', 'normal' or 'verbose' (do plot)
%
% OPTIONS OUT:
%  iter (integer): number of iterations processed.
%  kvg (boolean) : convergence flag (1 if ok).
%
% Use 'leasqrexamp' for a fitting example.

% Author:  EF <manuf@ldv.univ-montp2.fr>, RM <rfm2@ds2.uh.cwru.edu>, AJ <jutan@charon.engga.uwo.ca>
% Description: Non Linear Least Square multivariable fit.

% --------------------------- Function Calls ------------------------
% dfdp.m   : function file computing partials by finite differencies. Do not change.
% leasqrfunc.m   : example of function f(x,p) that can be used.
% leasqrexamp.m  : fit exemple with Leasqr V3.b

% Adapted for Octave : E.Farhi.   04/98   (manuf@ldv.univ-montp2.fr)
% Richard I. Shrager (301)-496-1122
% Modified by A.Jutan (519)-679-2111  jutan@charon.engga.uwo.ca
% Modified by Ray Muzic 14-Jul-1992   rfm2@ds2.uh.cwru.edu 

% Version 3.beta    revised 07/97.
% Part of 'Spectral tools'. E.Farhi. 06/97

% Refrences:
% Bard, Nonlinear Parameter Estimation, Academic Press, 1974.
% Draper and Smith, Applied Regression Analysis, John Wiley and Sons, 1981.

% argument processing -----------------------------------------------

if ~exist('last')
	last = 'normal';
end

if  isempty(last)
	last = 'normal';
end

if strcmp(last,'verbose')
	tmp = 3;
elseif strcmp(last,'silent')
	tmp = 0;
elseif strcmp(last,'min')
	tmp = 1;
else
	tmp = 2;
end

if tmp>0
	fprintf(1,'\nMarquardt-Levenberg Multivariable Fit\n');
end

if ~exist('ldfdp') ldfdp=[]; end
if ~exist('dp') dp=[]; end
if ~exist('wt') wt=[]; end
if ~exist('niter') niter=[]; end
if ~exist('stol') stol=[]; end

if isempty(ldfdp)
	ldfdp='dfdp';
	if (tmp>1)
		disp('using file dfdp.m for partials computation');
	end
end						% partials
if isempty(dp)
	dp=.01;
	if (tmp>1)
		disp('Increments for partials are dp = 0.01*pin');
	end
end

if (length(dp) == 1)
	dp=dp*pin;
	n = find(dp == 0);
	if (~isempty(n))
		dp(n) = 0.001;
	end
end	
					% step
if isempty(wt)
	wt=1./sqrt(abs(ys+(ys == 0)));
	if (tmp>1)
		disp('using standard normalisation (1/sqrt(y))')
	end
end						% data weights
if isempty(niter)
	niter=50;
	if (tmp>1)
		fprintf(1,'Iterations : %i\n',niter);
	end
end						% max nb of iters
if isempty(stol)
	stol=.001;
	if (tmp>1)
		fprintf(1,' Tolerancy : %g\n',stol);
	end
end						% tolerancy

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
x=vect2column(x);
y=vect2column(ys);
wt=vect2column(wt);
pin=vect2column(pin);
dp=vect2column(dp);
% check data vectors- same length?
m=length(y); 
n=length(pin);
p=pin;[m1,m2]=size(x);
if (m1~=m)
	error('input(x)/output(y) data must have same number of rows ');
end

if ~exist('options') 
	options = [];
end
if isempty(options)
  options=[zeros(n,1) Inf*ones(n,1)];
  nor = n; 
  noc = 2;
else
  [nor noc]=size(options);
  if (nor ~= n),
    error('options and parameter matrices must have same number of rows'),
  end;
  if (noc ~= 2),
    options=[options(noc,1) Inf*ones(noc,1)];
  end
end
pprec=options(:,1);
maxstep=options(:,2);

if (m > 10000)
    warning('HUH ! That makes a lot of points ! I may crash during fit...');
end

if (m <= n)
        disp('Warning : need to extend data...');
        newm = ceil(n/m+1)*m;
        newx = linspace(x(1), x(length(x)), newm);
        newx = vect2column(newx);
        y = interpsp(x,y,newx);
        wt = interpsp(x,wt,newx);
        x=newx;
        m = newm;
end

% set up for iterations (with initial data) -------------------------

f=vect2column(feval(func,x,p)); 
fbest=f; pbest=p;		% initial function values
r=wt.*(y-f);					% error
sbest=r'*r;					% square error
nrm=zeros(n,1);
chgprev=Inf*ones(n,1);
kvg=0;
epsLlast=1;
epstab=[.1 1 1e2 1e4 1e6];

t0 = clock;

% do iterations -----------------------------------------------------

if tmp>0
	disp('')
	disp(sprintf('*Beginning fit (max %d iterations)',niter));
	if tmp>1
	disp('--------------------------------------')
	disp('Iteration  Time(s)  Residual Variance');
	end
end

tic
for iter=1:niter,
  pprev=pbest;					% current params
  prt=vect2column(feval(ldfdp,x,fbest,pprev,dp,func));	% partials of func for x.
  r=wt.*(y-fbest);				% error
  sprev=sbest;					% square error
  sgoal=(1-stol)*sprev;				% sq error to achieve
  for j=1:n,			% computes partials norm
    if dp(j)==0,
      nrm(j)=0;
    else
      prt(:,j)=wt.*prt(:,j);
      nrm(j)=prt(:,j)'*prt(:,j);
      if nrm(j)>0,
        nrm(j)=1/sqrt(nrm(j));
      end
    end
    prt(:,j)=nrm(j)*prt(:,j);	% normalizes partials
  end
  [prt,s,v]=svd(prt,0);		% prt=unit matr, s=eigenval : prt=prt*s*v'.
  s=diag(s);			% get eigenvalues (diagonal of sigma)
  g=prt'*r;			% gradient by Gauss-Newton formula.
  for jjj=1:length(epstab),
    epsL = max(epsLlast*epstab(jjj),1e-7);
    se=sqrt((s.*s)+epsL);
    gse=g./se;
    chg=((v*gse).*nrm);		% change on params
				% check the change constraints and apply as necessary
    ochg=chg;
    for iii=1:n,
      if (maxstep(iii)==Inf), break; end;
      chg(iii)=max(chg(iii),-abs(maxstep(iii)*pprev(iii)));
      chg(iii)=min(chg(iii),abs(maxstep(iii)*pprev(iii)));
    end;
    if (any(ochg ~= chg)),
      if (tmp>0)
        disp(['Change in parameter(s): ' ...
         sprintf(1,'%d ',find(ochg ~= chg)) 'were constrained']);
      end
    end;
    aprec=abs(pprec.*pbest);       %---
    if (any(abs(chg) > 0.1*aprec)) %---  % only worth evaluating function if
      p=chg+pprev;                        % there is some non-miniscule change
      f=vect2column(feval(func,x,p));
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
    if (tmp>0)
      fprintf(1,'Parameter changes converged to specified precision\n');
    end
    break;
  else
    chgprev=chg;
  end;
  if ss>sgoal,
    break;
  end;
  if tmp>1
    fprintf(1,'   %3d ',iter);
    fprintf(1,'      %6.2f   %8.3f\n', toc, ss/length(y));
  end

end;

% set return values -------------------------------------------------
%
p=pbest;
f=fbest;
ss=sbest;
kvg=((sbest>sgoal)|(sbest<=eps)|kvg);
if (tmp>0)
	if (kvg ~= 1)
		fprintf(1,'\n\n** CONVERGENCE NOT ACHIEVED! **  ');
	else
		fprintf(1,'\n\n-- The fit process has converged --  ');
	end;
	fprintf(1,'%3d iterations, %5.1f s.\n\n',iter,etime(clock, t0));
end

%Calculate R (Ref Draper & Smith p.46)
%
r=corrcoef(y.*wt,f.*wt);
r2=r.*r';
ss=sum(((f-y).*wt).^2)/length(y);

if (m > 500)	% try to avoid crash if too many points
  if (tmp>1)
  	disp('Reducing data for statistics to 500 points');
  end
  xsav = x;
  ysav = y;
  fsav = f;
  wsav = wt;
  k = m/500;
  k = floor(k*(1:500));
  x=x(k);
  y=y(k);
  f=f(k);
  if max(size(wt) > 1)
    wt = wt(k);
  end
  msav = m;
  m = 500;
end % if m

  % CALC VARIANCE COV MATRIX AND CORRELATION MATRIX OF PARAMETERS
  % re-evaluate the Jacobian at optimal values
  jac=vect2column(feval(ldfdp,x,f,p,dp,func));
  msk = find(dp ~= 0);
  n = length(msk);           	% reduce n to equal number of estimated parameters
  jac = jac(:, msk);		% use only fitted parameters

  %% following section is Ray Muzic's estimate for covariance and correlation
  %% assuming covariance of data is a diagonal matrix proportional to
  %% diag(1/wt.^2).  
  %% cov matrix of data est. from Bard Eq. 7-5-13, and Row 1 Table 5.1 

  Qinv=diag(wt.*wt);
  Q=diag((0*wt+1)./(wt.^2));

  resid=y-f;                                      %un-weighted residuals
%  Z=((m-n)*jac'*Qinv*jac)/(n*resid'*Qinv*resid); % Z= matrix that defines confidence region
  covr=resid'*Qinv*resid*Q/(m-n);                 %covariance of residuals
  Vy=1/(1-n/m)*covr;  % Eq. 7-13-22, Bard         %covariance of the data 
%  stdresid=resid./sqrt(diag(Vy));  % compute then convert for compact storage

  jtgjinv=pinv(jac'*Qinv*jac);
  if ~all(isnan(jtgjinv))
	covp=jtgjinv*jac'*Qinv*Vy*Qinv*jac*jtgjinv; % Eq. 7-5-13, Bard %cov of parm est
	corp = ones(n,n);
	for k=1:n,
		for j=k:n,
			if covp(k,k)*covp(j,j)
				corp(k,j)=covp(k,j)/sqrt(abs(covp(k,k)*covp(j,j)));
				corp(j,k)=corp(k,j);
			end
		end;
	end;
  else
    disp('Warn : leasqr : couldn''t compute squared Jacobian inverse matrix (no covp, no corp).');
    covp = zeros(n,n);
    corp = ones(n,n);
  end

if (exist('msav'))	% restoring data
  x=xsav;
  y=ysav;
  f=fsav;
  m=msav;
  wt = wsav;
end

deltap= sqrt(diag(covp));

j=1;
for i = 1:length(pin)
if (dp(i) == 0)
	pout(i) = pin(i);
	deltapout(i) = 0;
else
	pout(i) = p(i);
	deltapout(i) = deltap(j);
	j = j+1;
end
end % for i

%-------- Work out Chi squared ------------------------
% NB number of free pars isn't right: should be sum(1-dp) not length(p)
v=length(y)-length(find(dp));		% # of degrees of freedom (#sel points - #free pars)
ChiSq = sum( ((y-f).*wt).^2 )/v;

pout = vect2column(pout);
deltapout = vect2column(deltapout);

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
	if (tmp>1)
	disp('* Covariance matrix of estimated parameters')
	disp(covp)
	disp('* Correlation matrix of estimated parameters ')
	disp(corp)
	end
end

if (tmp>2)
  fprintf(1,'Fit results sent to graph...\n');
  plot(x(:,1),y,x(:,1),f);	% plot command
  xlabel('X'); ylabel('Y');
  title('Fit result');
end

if (size(ys,1)==1)
	f = f';
end

p=p';
