close all
clear all
clear gobal
%
% TRIXFIT script
%
% Fitting triple-axis data, including resolution effects calculated using the
% Cooper-Nathans or Popovici methods.
%
% Notes: 
% (1)  The Monte Carlo integration routine used for the convolution introduces
%      noise, and as a result it is necessary to balance the parameters used
%      by the fit routine, the three fcp parameters:
%                [parameter_step no_of_iterations convergence_criterion]
%      against the number of sampling points used by the Monte Carlo. In general,
%      the larger the number of sampling points the smaller the paramter_step
%      can be between iterations. Typical values to get started are 1000 and 0.05.
% (2)  This version of TRIXIFT can cope with scans of H, K, L or E, but not mixed scans.
%      In the event of scans involving multiple reciprocal lattice coordinates, edit 
%      the file 'trixfit.m'
% (3)  Before fitting EACH AND EVERY data set it is necessary to call trixfit_ini 
%
% Des McMorrow and Henrik Roennow
% Version: May 2003
%

initialise.monitor_flag=1;                   % 1 for count on monitor, 0 for count on time
initialise.monte_carlo_samples=500000;        % Monte Carlo steps for convolution integration    

initialise.resolution_method='rc_popma';     % rc_popma (Popovici) or rc_cnmat (Cooper-Nathans)
%initialise.rescal_pars='U:\Data_NoBackup\neutrons\in14\dpip_in14_2006\dpip_in14_2006.par';      % parameters for Cooper-Nathans (mandatory)
%initialise.popovici_pars='U:\Data_NoBackup\neutrons\in14\dpip_in14_2006\dpip_in14_2006.cfg';    % parameters for Popovici (optional)

initialise.rescal_pars='/afs/psi.ch/user/t/thielemann/Data_NoBackup/neutrons/in14/dpip_in14_2006/dpip_in14_2006.par';      % parameters for Cooper-Nathans (mandatory)
initialise.popovici_pars='/afs/psi.ch/user/t/thielemann/Data_NoBackup/neutrons/in14/dpip_in14_2006/dpip_in14_2006.cfg';    % parameters for Popovici (optional)

initialise.xsec_file='dpip_xsec_notest';            % definition of S(Q,w)
initialise.bkgd_file='dpip_bkgd';            % background definition 
initialise.pnam_file='dpip_pnam';            % sets parameter names
initialise.corr_file='lsco_corr';            % correction to calculated intensities, 
                                             % e.g. lambda/2 in monitor

% Define a few variables to make life easier

niter=1000;   % No. of iterations for fitting
pvals=[];  % Fitted parameters
evals=[];  % Errors in parameters
mon6d7min=6000; % mn cnts for dpip 2006 in14 experiment
DAT_PATH_IN14=char('/afs/psi.ch/user/t/thielemann/Data_NoBackup/neutrons/in14/dpip_in14_2006/data/');
Jl=0.284;

%======================================================================================                                             
                                             
%----- Run trifit_ini to initialise the parameters

error_status=trixfit_ini(initialise);
if ~isempty(error_status), disp('Error initialising parameters'),return, end
                    
%----- Load and combine data
s=loads('illbatch',[DAT_PATH_IN14 '042[370 383 384 420 421],X=QH,Y=CNTS,M=M1'])*mon6d7min;
s=combine(0.005,s);

np=length(getfield(s,'x'));
%pin=[1000 1.105 0 0.3]; %----- Specify Qh,Qk,Ql and w, 1000 indicates scan variable
pin=[0.0253 0 0.9947 0.2 1.025 0 0.185 0.2 np]; %----- Specify start and end points Qh1,Qk1,Ql1,w1 Qh2,Qk2,Ql2,w2 of scan
                 %----- and the number of points in the scan
pin=[pin 0.2 0.5 0 100 0.284];   %----- Cross-sec parameters defined in xsec_file
pin=[pin 200 0 0 0];            %----- Background parameters used in bkgd_file

np=100;
m=0.5;
pin2=[0.0 0 0.0 0.2 1 0 0 0.2 np]; %----- Specify start and end points Qh1,Qk1,Ql1,w1 Qh2,Qk2,Ql2,w2 of scan
                 %----- and the number of points in the scan
pin2=[pin2 0.2 m 0 1999.999 0.284];   %----- Cross-sec parameters defined in xsec_file
pin2=[pin2 200 0 0 0];            %----- Background parameters used in bkgd_file

%----- Set fit control parameters (fcp), determine which parameters to fit and perform fit  

fcp=[1 100 0.001];notfixed=zeros(size(pin));
%fcp=[0.01 niter 0.0001];notfixed=zeros(size(pin));
notfixed([13])=1;
[s,f]=fits(s,'trixfit',pin,notfixed,fcp)
np=40;
sigma=0.015;
xmin=min(getfield(s,'x'))
xmax=max(getfield(s,'x'));
x=[xmin:(xmax-xmin)/(np-1):xmax];

error_status=trixfit_ini(initialise);
if ~isempty(error_status), disp('Error initialising parameters'),return, end
y=trixfit(x,f.pvals);
%x=getfield(s,'x');
%y=trixfit(x,pin);
%y=trixfit(x,pin2);
A=x'*ones(1,length(x));
B=A';
C=exp(-(B-A).^2./sigma^2)./sqrt(2*pi*sigma);
y1=(y.*[diff(x) x(length(x))-x(length(x)-1)])*C
y1(1)=2.*y1(1);
y1(length(y1))=2.*y1(length(y1));
hold on
plot(x,y,'bo',x,y1/0.28,'ro')
plot(s)
%[s,f]=fits(s,'trixfit',pin,notfixed)
pvals=[pvals [f.pvals]]; evals=[evals [f.evals]];


% test if the fit routine works, with a simple cross section and a spec1d
% object that shall be fitted in:
np=50;
m=0.5;
pin2=[0.0253 0 0.9947 0.2 1.025 1 0.185 0.2 np]; %----- Specify start and end points Qh1,Qk1,Ql1,w1 Qh2,Qk2,Ql2,w2 of scan
                 %----- and the number of points in the scan
pin2=[pin2 0.2 m 0 9 0.284];   %----- Cross-sec parameters defined in xsec_file
pin2=[pin2 200 0 0 1];            %----- Background parameters used in bkgd_file
fcp=[0.05 100 0.000001];
notfixed=zeros(size(pin2));
notfixed([13])=1;
x=[-0.4:0.05:1.4];
sigma_w=0.1;
sigma_Q=0.1;
y_disp=0.2.*cos(2*pi.*x)+0.15;
y_qsc=exp(-((pin2(4)-0.1.*cos(2*pi.*x)-0.1)/sigma_w).^2-(x-1/(2*pi).*real(acos((y_disp-0.1)/0.1))).^2);
s=spec1d(x,y_qsc,sqrt(y_qsc)./10);

y=trixfit(x,pin2)
[s,f]=fits(s,'trixfit',pin2,notfixed,fcp)



  