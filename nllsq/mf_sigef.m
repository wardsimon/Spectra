function [sig,covp,corp,r2,rv,f] = mf_sig(x,y,wt,p,dp,notfixed,func,f)
% [sig,covp,corp,r2,rv,f] = mf_sig(x,y,wt,p,dp,notfixed,func,f)
% mf_sig computes sigma for each parameter
% notfixed is 0 for fixed parameter, 1 for not fixed
% good dp is 0.01

m=length(y);
x=x(:);y=y(:); wt=wt(:); p=p(:); notfixed = notfixed(:);
dp = dp(:); 

if (nargin <= 8) f=[]; end
if isempty(f) | (length(f) ~= m)
	f=feval(func,x,p);
end
f=f(:);
if length(wt) ~= length(y)
	wt = ones(length(y),1);
end
index = find(~isnan(x) & ~isnan(y) & ~isnan(wt)& ~isnan(f) & ~isinf(f) & ~isinf(y) & ~isinf(wt) & wt);
x=x(index);
y=y(index);
wt = wt(index);
f=f(index);
m=length(y);

dp = dp.*notfixed;


sig = []; covp = []; corp = []; r2 = []; rv = [];

if isempty(p)
	return
end

if (m > 500)	% try to avoid crash if too many points
  disp('Reducing data for statistics to 500 points');
  xsav = x;
  ysav = y;
  fsav = f;
  wsav = wt;
  k = m/500;
  k = floor(k*(1:500));
  k = k(find(k >= 1 & k <= m));

  x=x(k);
  y=y(k);
  f=f(k); 

  if (length(wt))
    wt = wt(k);
  end
  msav = m;
  m = 500;
end % if m

  % CALC VARIANCE COV MATRIX AND CORRELATION MATRIX OF PARAMETERS
  % re-evaluate the Jacobian at optimal values
  jac=feval('mf_dfdp',x,f,p,dp,func);
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
%  covr=resid'*Qinv*resid*Q/(m-n);                 %covariance of residuals
%  Vy=1/(1-n/m)*covr;  % Eq. 7-13-22, Bard         %covariance of the data 
%  stdresid=resid./sqrt(diag(Vy));  % compute then convert for compact storage

  jtgjinv=pinv(jac'*Qinv*jac);
  if (exist('jtgjinv'))
    covp=resid'*Qinv*resid/(m-n)*jtgjinv;
    d=sqrt(abs(diag(covp)));
    corp=covp./(d*d');
  else
    disp('Warn : couldn''t compute squared Jacobian inverse matrix (no covp, no corp).');
    covp = zeros(n,n);
    corp = ones(n,n);
  end


stdp=sqrt(diag(covp));

if (exist('msav'))	% restoring data
  x=xsav;
  y=ysav;
  f=fsav;
  m=msav;
  wt = wsav;
end

j=1;
sig=zeros(size(p));
for i=1:length(stdp)
	while notfixed(j)==0
		j=j+1;
	end
	sig(j)=stdp(i);
	j=j+1;
end

r=corrcoef(y.*wt,f.*wt);
r2=r.*r';
rv=sum(((f-y).*wt).^2/length(y));
