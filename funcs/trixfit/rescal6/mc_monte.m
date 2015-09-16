% This function calculates a Monte Carlo
% distribution of points in (q,w) space
% which is Gaussian w.r.t. the priciple 
% axes of the resolution ellipsoid.

function [X]=rc_monte(M,N)

% M :== the resolution ellipsoid in the hkl axes frame.
% N :== the number of random points required.

c=2.3548/(2*sqrt(2*log(2)));  % to convert to sigma.

[V,E]=eig(M);
sigma=zeros(1,4);
sigma(1)=c/sqrt(E(1,1));
sigma(2)=c/sqrt(E(2,2));
sigma(3)=c/sqrt(E(3,3));
sigma(4)=c/sqrt(E(4,4));

% Random number generator produces a Gaussian distribution with
% a sigma=1. The generated distribution is scaled by the sigma of
% the resolution ellipsoid principle axes.

xp(1,:)=sigma(1)*randn(1,N);
xp(2,:)=sigma(2)*randn(1,N);
xp(3,:)=sigma(3)*randn(1,N);
xp(4,:)=sigma(4)*randn(1,N);
X=inv((V'))*(xp);
