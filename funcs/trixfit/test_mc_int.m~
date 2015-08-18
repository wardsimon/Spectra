% m-file to check the monte-carlo integration
close all
clear all
clear gobal


mon6d7min=6000; % mn cnts for dpip 2006 in14 experiment
DAT_PATH_IN14=char('/afs/psi.ch/user/t/thielemann/Data_NoBackup/neutrons/in14/dpip_in14_2006/data/');
Jl=0.284;

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
                                             
                    
%----- Load and combine data
s=loads('illbatch',[DAT_PATH_IN14 '042[370 383 384 420 421],X=QH,Y=CNTS,M=M1'])*mon6d7min;
s=combine(0.005,s);                                            

np=length(getfield(s,'x'));
%pin=[1000 1.105 0 0.3]; %----- Specify Qh,Qk,Ql and w, 1000 indicates scan variable
pin=[0.0253 0 0.9947 0.58 1.025 0 0.185 0.58 np]; %----- Specify start and end points Qh1,Qk1,Ql1,w1 Qh2,Qk2,Ql2,w2 of scan
                 %----- and the number of points in the scan
pin=[pin 0.2 0.5 0 1e3 0.284];   %----- Cross-sec parameters defined in xsec_file
pin=[pin 200 0 0 0];            %----- Background parameters used in bkgd_file

xmin=min(getfield(s,'x'))
xmax=max(getfield(s,'x'));
x=[xmin:(xmax-xmin)/(np-1):xmax];

figure
hold on;
for j=1:1
    error_status=trixfit_ini(initialise);
    if ~isempty(error_status), disp('Error initialising parameters'),return, end
    y=trixfit(x,pin);
    plot(x,y,'Color',[j*0.20 j*0.1 j*0.1]);
  %  waitforbuttonpress
end

hold off;

initialise.monte_carlo_samples=50;        % Monte Carlo steps for convolution integration    
figure
hold on;
for j=1:10
    error_status=trixfit_ini(initialise);
    if ~isempty(error_status), disp('Error initialising parameters'),return, end
    y=trixfit(x,pin);
    plot(x,y,'Color',[j*0.10 j*0.07 j*0.05]);
end

hold off;

initialise.monte_carlo_samples=500;        % Monte Carlo steps for convolution integration    
figure
hold on;
for j=1:10
    error_status=trixfit_ini(initialise);
    if ~isempty(error_status), disp('Error initialising parameters'),return, end
    y=trixfit(x,pin);
    plot(x,y,'Color',[j*0.10 j*0.07 j*0.05]);
end

hold off;

initialise.monte_carlo_samples=5000;        % Monte Carlo steps for convolution integration    
figure
hold on;
for j=1:10
    error_status=trixfit_ini(initialise);
    if ~isempty(error_status), disp('Error initialising parameters'),return, end
    y=trixfit(x,pin);
    plot(x,y,'Color',[j*0.10 j*0.07 j*0.05]);
end

hold off;
initialise.monte_carlo_samples=50000;        % Monte Carlo steps for convolution integration    
figure
hold on;
for j=1:10
    error_status=trixfit_ini(initialise);
    if ~isempty(error_status), disp('Error initialising parameters'),return, end
    y=trixfit(x,pin);
    plot(x,y,'Color',[j*0.10 j*0.07 j*0.05]);
end

hold off;

initialise.monte_carlo_samples=1000000;        % Monte Carlo steps for convolution integration    
figure
hold on;
for j=1:10
    error_status=trixfit_ini(initialise);
    if ~isempty(error_status), disp('Error initialising parameters'),return, end
    y=trixfit(x,pin);
    plot(x,y,'Color',[j*0.10 j*0.07 j*0.05]);
end

hold off;
