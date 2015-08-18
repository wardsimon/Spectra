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
%
% Output variables: s=calculated intensity
%


%----- Edit Cross-section definition after this point
%global grid_caux;
global spm;
global smp;
global szz;
global k;
global e;

%DATA_PATH_NUM_CAUX_NEW='U:\numerics\caux\';
%DATA_PATH_NUM_CAUX_NEW='/afs/psi.ch/user/t/thielemann/numerics/caux/';
DATA_PATH_NUM_CAUX_NEW='D:\PhD\data\numerics\caux\field_scan\';
DATA_PATH_NUM_CAUX='D:\PhD\data\numerics\caux\';

m=p(11); % m is the relative ladder magnetisation
dsf=p(12);
Jl=p(14);
a1=[0.3904;-0.1598;0.4842];
a2=[0.3904;0.1598;0.4842];

ar=[0.76065608444616;0;0]; % the reciprocal lattice vectors
br=[0;0.37003447038749;0];
cr=[0.08301979445134;0;0.50752708458640];

A_0=0.0232; % the constants needed for the copper magnetic form factor. It will be calculated with
a_0=34.9686;% this m-file since the function call of cu_mff costs a lot of time compared to the computation
B_0=0.4023;% it does.
b_0=11.5640;
C_0=0.5882;
c_0=3.8428;
D_0=-0.0137;


if load_the_matrices==1
    if m < 0.5
        k=load([DATA_PATH_NUM_CAUX_NEW 'Kfull_N_200.dat']);
        e=load([DATA_PATH_NUM_CAUX_NEW 'omega_Nw_1024_wmax_5.dat']);
    
        smp=load([DATA_PATH_NUM_CAUX_NEW 'Smp_JSC_D_0.500001_N_200_M_' num2str(m*200) '_Nw_1024_wmax_5_w_0.02_bin_sm.dat']);
        spm=load([DATA_PATH_NUM_CAUX_NEW 'Spm_JSC_D_0.500001_N_200_M_' num2str(m*200) '_Nw_1024_wmax_5_w_0.02_bin_sm.dat']);
        szz=load([DATA_PATH_NUM_CAUX_NEW 'Szz_JSC_D_0.500001_N_200_M_' num2str(m*200) '_Nw_1024_wmax_5_w_0.02_bin_sm.dat']);
%        grid_caux=s_long+0.25.*smp+0.25.*spm;
        disp(['loaded m = ' num2str(m) ' cross-section']);
    elseif m == 0.5
         k=load([DATA_PATH_NUM_CAUX 'K_N_200.dat']);
         e=load([DATA_PATH_NUM_CAUX 'omega_N_200.dat']);
         smp=load([DATA_PATH_NUM_CAUX 'Smp_JS_D_0.4999_N_200_M_100_hmax_3.dat']);
         spm=smp;
         szz=load([DATA_PATH_NUM_CAUX 'Szz_JS_D_0.4999_N_200_M_100_hmax_3.dat']);
%         grid_caux=(s_long+s_trans);
    end
end

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
% restrict the Caux-grid to reasonable size, such that less operations on
% zeros are performed:

e=e(1:find(e==e(3.998<e & e<3.998+e(2)-e(1))));
kk=k/2/pi;
%Qhh=abs(mod(Qh,1));
%Qhh=mod(abs(Qh),1);
Qhh=mod(0.981.*Qh,1);
iQhh=Qhh<0;
Qhh(iQhh)=Qhh(iQhh)+1;
ik=ones(size(Qh));
%dik=inf(size(Qh));

dk=diff(kk);
dk=dk(1);

ik=ceil(abs(Qhh./dk));
  
ee=e.*Jl;
ie=ones(size(w));
de=diff(ee);
de=de(1);

ie=ceil(abs(w./de));

% n_w(find(ie~=1))=1./(exp(w(find(ie~=1)).*(11.609/p(10)))-1);
% n_w(find(ie==1))=0;

s=0;
gclen=length(szz(:,1));
%gclen=length(grid_caux(:,1));
if dsf==0
    for j=1:length(n_w)
        if ie(j) < 4
            s=s+1;
        else
%            s=s+(n_w(j)+1).*sign(w(j)).*(grid_caux(ie(j),ik(j)));
            s=s+1;
        end
    end
else
%   figure
%     clf
%     hold on
     for j=1:length(n_w)
        if (ie(j) < 2) | (ie(j) > gclen)
   %         disp(['E-cut off']);
        else
            h_=Qh(j);
            k_=Qk(j);
            l_=Ql(j);
            Q_=sqrt((h_*ar(1)+l_*cr(1)).^2+(h_*ar(2)+l_*cr(2)).^2+(h_*ar(3)+l_*cr(3)).^2)/4/pi;
%            Q_=norm(h_.*ar+k_.*br+l_.*cr)/4/pi;
            mff=A_0.*exp(-a_0.*Q_.^2)+B_0.*exp(-b_0.*Q_.^2)+C_0.*exp(-c_0.*Q_.^2)+D_0;
%            mff=1;

%%%            s=s+(mff.^2).*(n_w(j)+1).*sign(w(j)).*(1-0.5*cos(2*pi*(h_*a1(1)+k_*a1(2)+l_*a1(3)))-0.5*cos(2*pi*(h_*a2(1)+k_*a2(2)+l_*a2(3))))*grid_caux(ie(j),ik(j))./w(j);
%            s=s+(mff.^2).*(n_w(j)+1).*sign(w(j)).*(1-0.5*cos(2*pi*(h_*a1(1)+k_*a1(2)+l_*a1(3)))-0.5*cos(2*pi*(h_*a2(1)+k_*a2(2)+l_*a2(3))))*grid_caux(ie(j),ik(j));
            s=s+(mff.^2).*(n_w(j)+1).*sign(w(j)).*((0.5-0.25*cos(2*pi*(h_*a1(1)+k_*a1(2)+l_*a1(3)))-0.25*cos(2*pi*(h_*a2(1)+k_*a2(2)+l_*a2(3))))*(spm(ie(j),ik(j))+smp(ie(j),ik(j)))...
              +(1+0.5*cos(2*pi*(h_*a1(1)+k_*a1(2)+l_*a1(3)))+0.5*cos(2*pi*(h_*a2(1)+k_*a2(2)+l_*a2(3)))).*szz(ie(j),ik(j)));
%            plot(h_,l_,'ro')
        end
    end
%      waitforbuttonpress
%      hold off
    s=s*0.5;
%    s=s*0.5;
end
 
s=1e7*p(13)./NMC.*s;% use this prefactor when working without Zheludev-det(M) trixfit.m
%====================================================================
function [Qh,Qk,Ql,w]=trixfit_fullmonte(NMC,Qvec,b_mat,fwhm)
%
% TRIXFIT function to generate Monte Carlo sampling of resolution volume
%
% Des McMorrow May 2003


xp=zeros(4,NMC);
%randn('state',0)
xp(1,:)=fwhm(1)*randn(1,NMC);
xp(2,:)=fwhm(2)*randn(1,NMC);
xp(3,:)=fwhm(3)*randn(1,NMC);
xp(4,:)=fwhm(4)*randn(1,NMC);
XMC=reshape(b_mat(1:16),4,4)'*xp;
 
Qh=XMC(1,:)+Qvec(1);% 1 x NMC matrix
Qk=XMC(2,:)+Qvec(2);% 1 x NMC matrix
Ql=XMC(3,:)+Qvec(3);% 1 x NMC matrix
w=XMC(4,:)+Qvec(4);% 1 x NMC matrix
