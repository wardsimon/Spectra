function [s]=lsco_xsec(NMC,Qvec,p,fwhm,b_mat)
%
% TRIXFIT function to define S(Q,w) and perform Monte Carlo integral
%         over resolution function
%
% Des Mcmorrow, October 2001
%
% Units: At the moment will only work with meV
%
% Input variables: 
%
%        p(1:4) = Reserved!
%        p(5)   = Temperature (K)
%        p(6)   = qh zone centre 
%        p(7)   = qk qzone centre 
%        p(8)   = Intensity scale factor
%        p(9)   = Damping in Q
%        p(10)  = Incommensurate wavevector 
%
% Output variables: s=calculated intensity
%

%----- Generate Monte Carlo sample of resolution volume

[Qh,Qk,Ql,w]=trixfit_fullmonte(NMC,Qvec,b_mat,fwhm);

%----- Edit Cross-section definition after this point

%----- Temperature factor

n_w=1./(exp(w*11.609/p(5))-1);

qh=(Qh-p(6)); 
qk=(Qk-p(7)); 

delta=[p(10) p(10)];
s_qw1=p(9)^4./(qh.^2+qk.^2-2*qh*delta(1)-2*qk*delta(2)+delta(1)^2+delta(2)^2+p(9)^2).^2;

delta=[-p(10) p(10)];
s_qw2=p(9)^4./(qh.^2+qk.^2-2*qh*delta(1)-2*qk*delta(2)+delta(1)^2+delta(2)^2+p(9)^2).^2;

delta=[p(10) -p(10)];
s_qw3=p(9)^4./(qh.^2+qk.^2-2*qh*delta(1)-2*qk*delta(2)+delta(1)^2+delta(2)^2+p(9)^2).^2;

delta=[-p(10) -p(10)];
s_qw4=p(9)^4./(qh.^2+qk.^2-2*qh*delta(1)-2*qk*delta(2)+delta(1)^2+delta(2)^2+p(9)^2).^2;

s=sum(p(8)*(s_qw1+s_qw2+s_qw3+s_qw4).*(n_w+1))/NMC;

s=s*2.9e13;

%====================================================================
function [Qh,Qk,Ql,w]=trixfit_fullmonte(NMC,Qvec,b_mat,fwhm)
%
% TRIXFIT function to generate Monte Carlo sampling of resolution volume
%
% Des McMorrow May 2003

xp=zeros(4,NMC);

xp(1,:)=fwhm(1)*randn(1,NMC);
xp(2,:)=fwhm(2)*randn(1,NMC);
xp(3,:)=fwhm(3)*randn(1,NMC);
xp(4,:)=fwhm(4)*randn(1,NMC);
XMC=reshape(b_mat(1:16),4,4)'*xp;

Qh=XMC(1,:)+Qvec(1);
Qk=XMC(2,:)+Qvec(2);
Ql=XMC(3,:)+Qvec(3);
w=XMC(4,:)+Qvec(4);