function [rho]=rc_focus(p)
%
% MATLAB function to calculate the optimum setting of the
% curvatures for the monochromator and analyser.
%
% input : parameters read from rescal windows
% output: rho returns the matrix of curvatures
%
% DFM 14.5.96
% DFM 7.8.2002
% Simon Ward 7.8.2013 Looked a bit wrong! 

verbose=0;
f=0.4826;

%----- Download rescal parameters from windows

if nargin ==0

p=rc_savp('respars');
verbose=1;

end

%----- Spectrometer distances
% These were cm and are converted to m
L1=p(62)/100; % Length Source->Mono
L2=p(63)/100; % Length Mono->Sample
L3=p(64)/100; % Length Sample->Analyser
L4=p(65)/100; % Length Analyser->Detector

%----- Angles at monochromator and analyser

dm=p(1);            % monochromator d-spacing in Angs.
da=p(2);            % analyser d-spacing in Angs.
kfix=p(9);          % fixed momentum component in ang-1.
fx=p(10);           % fx=1 for fixed incident and 2 for scattered wavevector.
w=p(34);            % energy transfer.

ki=sqrt(kfix^2+(fx-1)*f*w);  % kinematical equations.
kf=sqrt(kfix^2-(2-fx)*f*w);

theta_a=asin(pi/(da*kf));      % theta angles for analyser
theta_m=asin(pi/(dm*ki));      % and monochromator.

L= 1/(1/L1+1/L2);
rho_mv=abs(2*L*sin(theta_m));
rho_mh=abs(2*L/sin(theta_m));
L = 1/(1/L3+1/L4);
rho_av=abs(2*L*sin(theta_a));
rho_ah=abs(2*L/sin(theta_a));


% rho_mh=sin(theta_m)/(2*l_ms);
% rho_mv=1./(2*l_ms*sin(theta_m));
% rho_ah=(l_sa*sin(theta_a)+l_ad*sin(theta_a))/(2*l_sa*l_ad);
% rho_av=(1./l_sa+1./l_ad)/(2*sin(theta_a));

% The curvatures are in m^-1
rho=[rho_mh rho_mv rho_ah rho_av];

if verbose==1
 outstr=['Optimum focussing settings: [ ' num2str(rho)  ' ]'];
 disp(outstr);
 helpdlg(outstr,'Focussing settings');
end