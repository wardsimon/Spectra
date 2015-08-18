function [pnames]=trixfit_pnames(p)
%
% TRIXFIT function to define the fitting parameter names
%
% Des McMorrow, June 2001
%

pnames=strvcat('QH1','QK1','QL1','EN1','QH2','QK2','QL2','EN2','NP');
pnames=strvcat(pnames,'T(K)','rel. mag.','dsf on/off','Amp','Jl');
pnames=strvcat(pnames,'Bckgd','Slope','AmpE','ExpE','scan var.');

if length(pnames)~=length(p)
   warning('Number of parameter names not equal to number of parameters')
end