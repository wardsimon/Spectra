% This function calculates a grid
% distribution of points in (q,w) space
% which is Gaussian w.r.t. the priciple 
% axes of the resolution ellipsoid.

function [r,X]=rc_grid(M,INTLIM)

%  :== the resolution ellipsoid in the hkl axes frame.
%  :== the number of random points required.

c=2.3548/(2*sqrt(2*log(2)));  % to convert to sigma.

[V,E]=eig(M);
sigma=zeros(1,4);
sigma(1)=c/sqrt(E(1,1));
sigma(2)=c/sqrt(E(2,2));
sigma(3)=c/sqrt(E(3,3));
sigma(4)=c/sqrt(E(4,4));


% make up grid
int2=INTLIM^2;
int4=INTLIM^4;
a=ones(INTLIM,1)*[1:1:INTLIM];
arg1=ones(int2,1);
arg2=reshape(a,int2,1);
arg3=reshape(a',int2,1);
x1=reshape(arg1*arg2',int4,1);
x2=reshape(arg1*arg3',int4,1);
x3=reshape(arg2*arg1',int4,1);
x4=reshape(arg3*arg1',int4,1);


arg=2./INTLIM;
pr=exp(-(arg*[1:1:INTLIM]-0.5).^2/2);
pr=pr/sum(pr);

r=pr(x1).*pr(x2).*pr(x3).*pr(x4);
xp(1,:)=arg*sigma(1)*(x1-0.5)';
xp(2,:)=arg*sigma(2)*(x2-0.5)';
xp(3,:)=arg*sigma(3)*(x3-0.5)';
xp(4,:)=arg*sigma(4)*(x4-0.5)';

perm=[ -1 1 1 1 ; 1 -1 1 1; -1 -1 1 1; 1 1 1 -1;...
 -1 1 1 -1; 1 -1 1 -1; -1 -1 1 -1];

xt=[xp [-xp(1,:);xp(2:4,:)] [xp(1,:);-xp(2,:);xp(3:4,:)]...
 [-xp(1:2,:);xp(3:4,:)] [xp(1:3,:);-xp(4,:)] [-xp(1,:);xp(2:3,:);-xp(4,:)]...
 [xp(1,:);-xp(2,:);xp(3,:);-xp(4,:)] [-xp(1:2,:);xp(3,:);-xp(4,:)]];
r=[r r r r r r r r]/int4;
X=inv((V'))*(xt);
