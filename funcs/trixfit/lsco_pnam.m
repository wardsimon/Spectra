function [pnames]=trixfit_pnames(p)
%
% TRIXFIT function to define the fitting parameter names
%
% Des McMorrow, June 2001
%

pnames=strvcat('QH','QK','QL','EN');
pnames=strvcat(pnames,'T(K)','zc_H','zc_K','Amp','Gamma','Delta');
pnames=strvcat(pnames,'Bckgd','Slope','AmpE','ExpE');

if length(pnames)~=length(p)
   warning('Number of parameter names not equal to number of parameters')
end