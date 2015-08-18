function [y]=maxwell(E,T)
% calculates the normalized Maxwellian 
% distribution at an energy E (meV) for a thermalized
% source at a temperature T (K).
%
% Calling syntax:
%
% function [y]=maxwell(E,T)
%
% Alvis 12th December 1997
kb=1/11.8; % convert K to meV
Norm=1/(kb*T)*exp(-1);
Norm=1/Norm;
y=Norm*(E/(kb*T)^2).*exp(-E/(kb*T));
return
end
