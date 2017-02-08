function [s]=dpip_xsec_notest(NMC,Qvec,p,fwhm,b_mat,load_the_matrices)
%
% TRIXFIT function to define S(Q,w) and perform Monte Carlo integral
%         over resolution function
%
% Des Mcmorrow, October 2001
% Beni Thielemann, October 2006
%
% Units: At the moment will only work with meV
%
% Input variables: 
%
%        p(1:9) = Reserved!
%        p(10)   = Temperature (K)
%        p(11)   = relative magnetisation
%        p(12)   = dimer structure factor on=1/off=0
%        p(13)   = Intensity scale factor
%        p(14)   = Jl, energy scale of the continuum
%        p(15)   = sigma of gaussian linewidth
% Output variables: s=calculated intensity
%


%----- Edit Cross-section definition after this point
DATA_PATH_NUM_CAUX_NEW='D:\PhD\data\numerics\caux\';
dsf=p(12);
Jl=p(14);
sigma=p(11);
a1=[0.3904;-0.1598;0.4842];
a2=[0.3904;0.1598;0.4842];

ar=[0.75157718985402;0;0]; % the reciprocal lattice vectors
br=[0;0.37003447038749;0];
cr=[0.08301979445134;0;0.50752708458640];

A_0=0.0232; % the constants needed for the copper magnetic form factor. It will be calculated with
a_0=34.9686;% this m-file since the function call of cu_mff costs a lot of time compared to the computation
B_0=0.4023;% it does.
b_0=11.5640;
C_0=0.5882;
c_0=3.8428;
D_0=-0.0137;


 % now calculate s for a given Qh:

%----- Generate Monte Carlo sample of resolution volume
%[Qh,Qk,Ql,w]=trixfit_fullmonte(NMC,Qvec,b_mat,fwhm);


xp=zeros(4,NMC);
%randn('state',0)
xp(1,:)=fwhm(1)*randn(1,NMC);% create the Monte Carlo points with the correct
xp(2,:)=fwhm(2)*randn(1,NMC);% standard deviation
xp(3,:)=fwhm(3)*randn(1,NMC);
xp(4,:)=fwhm(4)*randn(1,NMC);
XMC=reshape(b_mat(1:16),4,4)'*xp;
 
Qh=XMC(1,:)+Qvec(1);% 1 x NMC matrix
Qk=XMC(2,:)+Qvec(2);% 1 x NMC matrix
Ql=XMC(3,:)+Qvec(3);% 1 x NMC matrix
w=XMC(4,:)+Qvec(4);% 1 x NMC matrix
% Qh,Qk,Ql,w are now Monte Carlo sampling points with coordinates in the
% eigenbasis of M. By picking the distribution of the points exactly like
% the exponential function exp(-0.5 Q' M Q), one only has to integrate the
% cross-section itself (Importance-sampling Monte-Carlo integration)

%----- Temperature factor
%n_w=1./(exp(w.*(11.609/p(10)))-1);% 1 x NMC matrix
Ec=0.001*p(10)/11.609; % this is the lowest cut-off for which no artefact around Qh = 1 appears.

n_w=1./(exp(w./Ec)-1).*(abs(w)>Ec)+sign(w)./(exp(Ec.*11.609./p(10))-1).*(abs(w)<Ec);% 1 x NMC matrix


s=0;
E_cutoff=0.001*Jl;
E_lowbound=pi/2*Jl*abs(sin(2*pi*Qh*0.987));

if dsf==0
%    s=sum((n_w+1).*sign(w)./(sqrt(2*pi)*sigma).*exp(-0.5.*((w-E_lowbound)/sigma).^2)./E_lowbound);
    s=sum((n_w+1).*sign(w)./(sqrt(2*pi)*sigma).*exp(-0.5.*((w-E_lowbound)/sigma).^2));
else
%     Q_=sqrt((h_*ar(1)+l_*cr(1)).^2+(h_*ar(2)+l_*cr(2)).^2+(h_*ar(3)+l_*cr(3)).^2)/4/pi;% 1/4pi factor is needed for mff calc. below
%     h_=Qh(:);
%     k_=Qk(:);
%     l_=Ql(:);
%     Q_=sqrt((h_.*ar(1)+l_.*cr(1)).^2+(h_.*ar(2)+l_.*cr(2)).^2+(h_.*ar(3)+l_.*cr(3)).^2)/4/pi;% 1/4pi factor is needed for mff calc. below
%     mff=A_0.*exp(-a_0.*Q_.^2)+B_0.*exp(-b_0.*Q_.^2)+C_0.*exp(-c_0.*Q_.^2)+D_0;
% %    s=sum((mff'.^2).*(n_w+1).*sign(w).*(1-0.5*cos(2*pi*(h_'.*a1(1)+k_'.*a1(2)+l_'.*a1(3)))-0.5*cos(2*pi*(h_'.*a2(1)+k_'.*a2(2)+l_'.*a2(3))))./(sqrt(2*pi)*sigma).*exp(-0.5.*((w-E_lowbound)/sigma).^2)./E_lowbound);
%     s=sum((mff'.^2).*(n_w+1).*sign(w).*(1-0.5*cos(2*pi*(h_'.*a1(1)+k_'.*a1(2)+l_'.*a1(3)))-0.5*cos(2*pi*(h_'.*a2(1)+k_'.*a2(2)+l_'.*a2(3))))./(sqrt(2*pi)*sigma).*exp(-0.5.*((w-E_lowbound)/sigma).^2));
% %    s=s+(mff.^2).*(n_w(j)+1).*sign(w(j)).*(1-0.5*cos(2*pi*(h_*a1(1)+k_*a1(2)+l_*a1(3)))-0.5*cos(2*pi*(h_*a2(1)+k_*a2(2)+l_*a2(3))))*grid_caux(ie(j),ik(j));
    for j=1:length(n_w)
        if w(j)> E_cutoff
            h_=Qh(j);
            k_=Qk(j);
            l_=Ql(j);
            Q_=sqrt((h_.*ar(1)+l_.*cr(1)).^2+(h_.*ar(2)+l_.*cr(2)).^2+(h_.*ar(3)+l_.*cr(3)).^2)/4/pi;% 1/4pi factor is needed for mff calc. below
            mff=A_0.*exp(-a_0.*Q_.^2)+B_0.*exp(-b_0.*Q_.^2)+C_0.*exp(-c_0.*Q_.^2)+D_0;
            %s=s+(mff.^2).*(n_w(j)+1).*sign(w(j)).*(1-0.5*cos(2*pi*(h_.*a1(1)+k_.*a1(2)+l_.*a1(3)))-0.5*cos(2*pi*(h_.*a2(1)+k_.*a2(2)+l_.*a2(3))))./(sqrt(2*pi)*sigma).*exp(-0.5.*((w(j)-E_lowbound(j))/sigma).^2)./w(j);
            s=s+(mff.^2).*(n_w(j)+1).*sign(w(j)).*(0.5-0.25*cos(2*pi*(0.987*h_*a1(1)+k_*a1(2)+l_*a1(3)))-0.25*cos(2*pi*(0.987*h_*a2(1)+k_*a2(2)+l_*a2(3))))...
              ./(sqrt(2*pi)*sigma).*exp(-0.5.*((w(j)-E_lowbound(j))/sigma).^2)./w(j);
            
            %            s=s+(mff.^2).*(n_w(j)+1).*sign(w(j)).*(1-0.5*cos(2*pi*(h_.*a1(1)+k_.*a1(2)+l_.*a1(3)))-0.5*cos(2*pi*(h_.*a2(1)+k_.*a2(2)+l_.*a2(3))))./(sqrt(2*pi)*sigma).*exp(-0.5.*((w(j)-E_lowbound(j))/sigma).^2);
        end
    end
end
 
s=1e7*p(13)./NMC.*s;%
