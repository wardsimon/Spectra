function [sout]=fits(varargin)
%%
% How to use:
% [yfit pin Chisquared]=fits(x,y,e,function,pin,fixedPin);
% fixedPin, 0=fixed, 1=free
x=varargin{1};
y=varargin{2};
e=varargin{3};
func=varargin{4};
pin=varargin{5};
notfixed=varargin{6};
fcp=[0.0001 20 0.0001];
options=[];

%----- Remove zeros from e
ezeros=find(e==0);
x(ezeros)=[]; y(ezeros)=[]; e(ezeros)=[];
%----- Fit data
[yfit,p,cvg,iter,corp,covp,covr,stdresid,Z,ChiSq,sig]=speclsqr(x,y,e,pin,notfixed,func,fcp,options);
sout{1}=yfit;
sout{2}=p;
sout{3}=ChiSq;