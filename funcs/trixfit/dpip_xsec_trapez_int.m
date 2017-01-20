function [s]=dpip_xsec_notest(NMC,Qvec,p,fwhm,b_mat)
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
global grid_caux;
global k;
global e;

%DATA_PATH_NUM_CAUX_NEW='U:\numerics\caux\';
DATA_PATH_NUM_CAUX_NEW='/afs/psi.ch/user/t/thielemann/numerics/caux/';

m=p(11); % m is the relative ladder magnetisation
dsf=p(12);
Jl=p(14);
a1=[0.3904;-0.1598;0.4842];
a2=[0.3904;0.1598;0.4842];

% load the appropriate numerical grid, if this hasn't been done so far
% Attention: if dpip_xsec has been called and loaded a grid and you want to
% use the cross section with a different m, it won't load the data with the
% new magnetisation until you do a grid_caux=[] !????????????
% ?????????????????????
 if isempty(grid_caux)                        
     k=load([DATA_PATH_NUM_CAUX_NEW 'Kfull_N_320.dat']);
     e=load([DATA_PATH_NUM_CAUX_NEW 'omega_Nw_1024_wmax_5.dat']);
     if m == 0.125
         smp0d125=load([DATA_PATH_NUM_CAUX_NEW 'Smp_JSC_D_0.5_N_320_M_40_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
         spm0d125=load([DATA_PATH_NUM_CAUX_NEW 'Spm_JSC_D_0.5_N_320_M_40_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
         szz0d125=load([DATA_PATH_NUM_CAUX_NEW 'Szz_JSC_D_0.5_N_320_M_40_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
         grid_caux=szz0d125+0.25.*smp0d125+0.25.*spm0d125;
     end
     if m == 0.25
         smp0d25=load([DATA_PATH_NUM_CAUX_NEW 'Smp_JSC_D_0.5_N_320_M_80_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);        
         spm0d25=load([DATA_PATH_NUM_CAUX_NEW 'Spm_JSC_D_0.5_N_320_M_80_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
         szz0d25=load([DATA_PATH_NUM_CAUX_NEW 'Szz_JSC_D_0.5_N_320_M_80_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
         grid_caux=szz0d25+0.25.*smp0d25+0.25.*spm0d25;
     end
     if m == 0.375
         smp0d375=load([DATA_PATH_NUM_CAUX_NEW 'Smp_JSC_D_0.5_N_320_M_120_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
         spm0d375=load([DATA_PATH_NUM_CAUX_NEW 'Spm_JSC_D_0.5_N_320_M_120_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
         szz0d375=load([DATA_PATH_NUM_CAUX_NEW 'Szz_JSC_D_0.5_N_320_M_120_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
         grid_caux=szz0d375+0.25*(smp0d375+spm0d375);
     end
     if m == 0.5
         k=load([DATA_PATH_NUM_CAUX_NEW 'K_N_200.dat']);
         e=load([DATA_PATH_NUM_CAUX_NEW 'omega_N_200.dat']);
         smp0d5=load([DATA_PATH_NUM_CAUX_NEW 'Smp_JS_D_0.4999_N_200_M_100_hmax_3.dat']);
         szz0d5=load([DATA_PATH_NUM_CAUX_NEW 'Szz_JS_D_0.4999_N_200_M_100_hmax_3.dat']);
         grid_caux=(szz0d5+0.5*smp0d5);
     end
 end

 % now calculate s for a given Qh:

xp=zeros(4,NMC);
xp(1,:)=[-fwhm(1):(2*fwhm(1)/(NMC-1)):fwhm(1)];
xp(2,:)=[-fwhm(2):(2*fwhm(2)/(NMC-1)):fwhm(2)];
xp(3,:)=[-fwhm(3):(2*fwhm(3)/(NMC-1)):fwhm(3)];
xp(4,:)=[-fwhm(4):(2*fwhm(4)/(NMC-1)):fwhm(4)];
XMC=reshape(b_mat(1:16),4,4)'*xp;
 
Qh=XMC(1,:)+Qvec(1);% 1 x NMC matrix
Qk=XMC(2,:)+Qvec(2);% 1 x NMC matrix
Ql=XMC(3,:)+Qvec(3);% 1 x NMC matrix
w=XMC(4,:)+Qvec(4);% 1 x NMC matrix

%----- Temperature factor
n_w=1./(exp(w.*(11.609/p(10)))-1);

% restrict the Caux-grid to reasonable size, such that less operations on
% zeros are performed:
e=e(1:find(e==e(3.996<e & e<4.005)));
kk=k/2/pi;
Qhh=abs(mod(Qh,1));

ik=ones(size(Qh));
dk=kk(2)-kk(1);
ik=ceil(Qhh./dk);
   
ee=e.*Jl;
ie=ones(size(w));
de=ee(2)-ee(1);

w=abs(w);
ie=ceil(w./de);

de=w(2)-w(1);
dk=Qhh(1)-Qhh(2);
s=0;

if dsf==0
    for j=1:length(n_w)
        s=s+(n_w(j)+1).*(grid_caux(ie(j),ik(j)))*de*dk;
%  s=sum(exp(-((Qhh-0.5)/0.1).^2))/NMC;
%        s=1;
% these lines are for testing if the fit routine works:
%        waitforbuttonpress;
%        length(n_w)
%        s=sum(p(13)*(n_w(j)+1).*exp(-((p(4)-0.2.*cos(2*pi.*Qh)-0.15)/sigma_w).^2-(Qh-1/(2*pi).*real(acos((0.2.*cos(2*pi.*Qh)-0.15)/sigma_Q))).^2))
    end
else 
    s=p(13)*sum(((1-0.5*cos(2*pi*(Qh*a1(1)+Qk*a1(2)+Ql*a1(3)))-0.5*cos(2*pi*(Qh*a2(1)+Qk*a2(2)+Ql*a2(3))))*grid_caux(ie,ik)).*(n_w+1))/NMC;
end
 
s=p(13).*2.9e16.*s;

