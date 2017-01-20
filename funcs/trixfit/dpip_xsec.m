function [s]=dpip_xsec(NMC,Qvec,p,fwhm,b_mat)
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

% % load the appropriate numerical grid, if this hasn't been done so far
% % Attention: if dpip_xsec has been called and loaded a grid and you want to
% % use the cross section with a different m, it won't load the data with the
% % new magnetisation until you do a grid_caux=[] !
% if isempty(grid_caux)                        
%     k=load([DATA_PATH_NUM_CAUX_NEW 'Kfull_N_320.dat']);
%     e=load([DATA_PATH_NUM_CAUX_NEW 'omega_Nw_1024_wmax_5.dat']);
%     if m == 0.125
%         smp0d125=load([DATA_PATH_NUM_CAUX_NEW 'Smp_JSC_D_0.5_N_320_M_40_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
%         spm0d125=load([DATA_PATH_NUM_CAUX_NEW 'Spm_JSC_D_0.5_N_320_M_40_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
%         szz0d125=load([DATA_PATH_NUM_CAUX_NEW 'Szz_JSC_D_0.5_N_320_M_40_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
%         grid_caux=2.9e13.*szz0d125+0.25.*smp0d125+0.25.*spm0d125;
%     end
%     if m == 0.25
%         smp0d25=load([DATA_PATH_NUM_CAUX_NEW 'Smp_JSC_D_0.5_N_320_M_80_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);        
%         spm0d25=load([DATA_PATH_NUM_CAUX_NEW 'Spm_JSC_D_0.5_N_320_M_80_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
%         szz0d25=load([DATA_PATH_NUM_CAUX_NEW 'Szz_JSC_D_0.5_N_320_M_80_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
%         grid_caux=2.9e13.*szz0d25+0.25.*smp0d25+0.25.*spm0d25;
%     end
%     if m == 0.375
%         smp0d375=load([DATA_PATH_NUM_CAUX_NEW 'Smp_JSC_D_0.5_N_320_M_120_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
%         spm0d375=load([DATA_PATH_NUM_CAUX_NEW 'Spm_JSC_D_0.5_N_320_M_120_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
%         szz0d375=load([DATA_PATH_NUM_CAUX_NEW 'Szz_JSC_D_0.5_N_320_M_120_Nw_1024_wmax_5_w_0.018_bin_sm.dat']);
%         grid_caux=2.9e13.*szz0d375+1/4*(smp0d375+spm0d375);
%     end
%     if m == 0.5
%         k=load([DATA_PATH_NUM_CAUX_NEW 'K_N_200.dat']);
%         e=load([DATA_PATH_NUM_CAUX_NEW 'omega_N_200.dat']);
%         smp0d5=load([DATA_PATH_NUM_CAUX_NEW 'Smp_JS_D_0.4999_N_200_M_100_hmax_3.dat']);
%         szz0d5=load([DATA_PATH_NUM_CAUX_NEW 'Szz_JS_D_0.4999_N_200_M_100_hmax_3.dat']);
%         grid_caux=2.9e13.*(szz0d5+0.5*smp0d5);
%     end
% end

% now calculate s for a given Qh:

%----- Generate Monte Carlo sample of resolution volume
%[Qh,Qk,Ql,w]=trixfit_fullmonte(NMC,Qvec,b_mat,fwhm);
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

%----- Temperature factor
n_w=1./(exp(w.*(11.609/p(10)))-1);

% % restrict the Caux-grid to reasonable size, such that less operations on
% % zeros are performed:
% e=e(1:find(e==e(3.996<e & e<4.005)));
% kk=k/2/pi;
% Qhh=mod(Qh,1);
% ik=ones(size(Qh));
% dik=inf(size(Qh));
% 
% L=length(k);
% for j=1:L
%     dnew=abs(Qhh-kk(j));
%     iii=dnew<dik;
%     ik(iii)=j;
%     dik(iii)=dnew(iii);
% end
% 
% % ie=zeros(1,length(Qhh));
% % L=length(k);
% % for j=1:L-1
% %     ie=ie +j.*((Qhh>k(j)) & (Qhh<k(j+1)));
% % end
%     
% ee=e.*Jl;
% ie=ones(size(w));
% die=inf(size(w));
% 
% L=length(ee);
% for j=1:L
%     dnew=abs(w-ee(j));
%     iii=dnew<die;
%     ie(iii)=j;
%     die(iii)=dnew(iii);
% end

sigma_w=0.1;
sigma_Q=0.1;
    

if dsf==0
    for j=1:length(n_w)
%        s=sum(p(13)*(n_w(j)+1).*(grid_caux(ie(j),ik(j))))/NMC;
% these lines are for testing if the fit routine works:
%        waitforbuttonpress;
%        length(n_w)
        s=sum(p(13)*(n_w(j)+1).*exp(-((p(4)-0.2.*cos(2*pi.*Qh)-0.15)/sigma_w).^2-(Qh-1/(2*pi).*real(acos((0.2.*cos(2*pi.*Qh)-0.15)/sigma_Q))).^2))
    end
else 
    s=sum(p(13)*((1-0.5*cos(2*pi*(Qh*a1(1)+Qk*a1(2)+Ql*a1(3)))-0.5*cos(2*pi*(Qh*a2(1)+Qk*a2(2)+Ql*a2(3))))*grid_caux(ie,ik)).*(n_w+1))/NMC;
end
 

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
