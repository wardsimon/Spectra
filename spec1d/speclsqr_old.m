function [p, std]=speclsqr(x,y,err,pin,dpin,func,fcp)
%
% Version 3.beta
% Levenberg-Marquardt nonlinear regression of f(x,p) to y(x)
% Richard I. Shrager (301)-496-1122
% Modified by A.Jutan (519)-679-2111
% Modified by Ray Muzic 14-Jul-1992

ss=inf; % matlab 6.5 started complaining about ss not being defined by compilation.

quiet=1;

if ~exist('fcp') | isempty(fcp)
	fcp=[0.0001 20 0.001];
end

dp=-dpin*fcp(1); niter=fcp(2); stol=fcp(3);

% global wt

% if isempty(wt)    
wt=1./err;
    

% end

if ~quiet
 disp('')
 disp(sprintf('*Beginning fit (max %d iterations)',niter));
 disp('--------------------------------------')
 disp('Iteration  Time(s)  Reduced Chi^2');
end
tic

y=y(:); wt=wt(:); pin=pin(:); dp=dp(:);
m=length(y); n=length(pin);

options=[zeros(n,1) Inf*ones(n,1)];
nor = n; noc = 2;
pprec=options(:,1);
maxstep=options(:,2);
p=pin;

f=feval(func,x,p);
fbest=f;
pbest=p;
% size(y)
% size(f)
% size(wt)
r=wt.*(y-f);
sbest=r'*r;
nrm=zeros(n,1);
chgprev=Inf*ones(n,1);
kvg=0;
epsLlast=1;
epstab=[.1 1 1e2 1e4 1e6];

% do iterations
for iter=1:niter,
  pprev=pbest;
  prt=feval('specdfdp',x,fbest,pprev,dp,func);
  r=wt.*(y-fbest);
  sprev=sbest;
  sgoal=(1-stol)*sprev;
  for j=1:n,
    if dp(j)==0,
      nrm(j)=0;
    else
      prt(:,j)=wt.*prt(:,j);
      nrm(j)=prt(:,j)'*prt(:,j);
      if nrm(j)>0,
        nrm(j)=1/sqrt(nrm(j));
      end;
    end
    prt(:,j)=nrm(j)*prt(:,j);
  end;
  [prt,s,v]=svd(prt,0);
  s=diag(s);
  g=prt'*r;
  for jjj=1:length(epstab),
    epsL = max(epsLlast*epstab(jjj),1e-7);
    se=sqrt((s.*s)+epsL);
    gse=g./se;
    try
        chg=((v*gse).*nrm);
    catch
        error('Number of fit parameters > Number of points')
        return
    end
%   check the change constraints and apply as necessary
    ochg=chg;
    for iii=1:n,
      if (maxstep(iii)==Inf), break; end;
      chg(iii)=max(chg(iii),-abs(maxstep(iii)*pprev(iii)));
      chg(iii)=min(chg(iii),abs(maxstep(iii)*pprev(iii)));
    end;
 
    aprec=abs(pprec.*pbest);
    if (any(abs(chg) > 0.1*aprec)),
      p=chg+pprev;
      f=feval(func,x,p);
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
    end;
  end;
  epsLlast = epsL;
  if ~exist('ss')
     break;
  elseif ss<eps,
    break;
  end
  aprec=abs(pprec.*pbest);
%  [aprec chg chgprev]
  if (all(abs(chg) < aprec) & all(abs(chgprev) < aprec)),
    kvg=1;
    break;
  else
    chgprev=chg;
  end;
  if ss>sgoal,
      break;
  end;
if ~quiet  
  disp(sprintf('   %3d      %6.2f   %8.3f', iter, toc, ss/(length(x)-sum(dpin))));
  disp(p')
end  
end;

% set return values
%
p=pbest;
f=fbest;
ss=sbest;
kvg=((sbest>sgoal)|(sbest<=eps)|kvg);
if kvg ~= 1 , disp(' CONVERGENCE NOT ACHIEVED! '), end;

% CALC VARIANCE COV MATRIX AND CORRELATION MATRIX OF PARAMETERS
% re-evaluate the Jacobian at optimal values
jac=feval('specdfdp',x,f,p,dp,func);
msk = dp ~= 0;
n = sum(msk);        % reduce n to equal number of estimated parameters
jac = jac(:, msk);	% use only fitted parameters

%% following section is Ray Muzic's estimate for covariance and correlation
%% assuming covariance of data is a diagonal matrix proportional to
%% diag(1/wt.^2).  
%% cov matrix of data est. from Bard Eq. 7-5-13, and Row 1 Table 5.1 

Qinv=diag(wt.*wt);
Q=diag((0*wt+1)./(wt.^2));
%[nrw ncw]=size(wt);
%Q=ones(nrw,ncw)./wt; Q=diag(Q.*Q);
resid=y-f;                                    %un-weighted residuals
covr=resid'*Qinv*resid*Q/(m-n);                 %covariance of residuals
Vy=1/(1-n/m)*covr;  % Eq. 7-13-22, Bard         %covariance of the data 
covr=diag(covr);                                %for compact storage
Z=((m-n)*jac'*Qinv*jac)/(n*resid'*Qinv*resid);
stdresid=resid./sqrt(diag(Vy));
jtgjinv=pinv(jac'*Qinv*jac);
covp=jtgjinv*jac'*Qinv*Vy*Qinv*jac*jtgjinv; % Eq. 7-5-13, Bard %cov of parm est
for k=1:n,
  for j=k:n,
    corp(k,j)=covp(k,j)/sqrt(abs(covp(k,k)*covp(j,j)));
    corp(j,k)=corp(k,j);
  end;
end;

std=sqrt(diag(covp));

j=1;
sig=zeros(size(p));
for i=1:length(std)
	while dp(j)==0
		j=j+1;
	end
	sig(j)=std(i);
	j=j+1;
end
std=sig;

%%% alt. est. of cov. mat. of parm.:(Delforge, Circulation, 82:1494-1504, 1990
%%disp('Alternate estimate of cov. of param. est.')
%%acovp=resid'*Qinv*resid/(m-n)*jtgjinv


