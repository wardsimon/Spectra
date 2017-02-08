function [s]=mc_sqw(disper,qh,qk,ql,w)
%
% MATLAB  function to define S(Q,w) for a sinusoidal dispersion surface
%
% Des Mcmorrow 7.11.95
% 
% Notes: (1) The lineshape is assumed to be Lorentzian
%        (2) S(Q,w) is assumed to follow a 3D dispersion
%
% Called by: mc_conv
% Calls  to: 
%
% Units: At the moment will only work with meV
%
% Input variables:      disp = dispersion parameters
%			disper(1)  = qh zone-boundary energy
%			disper(2)  = qk      "          "
% 			disper(3)  = ql      "          "
% 			disper(4)  = zone-centre energy gap
% 			disper(5)  = qh zone-centre
% 			disper(6)  = qk zone-centre
% 			disper(7)  = ql zone-centre
%                       disper(8)  = qh phase factor (multiples of pi)
%                       disper(9)  = qk phase factor (multiples of pi)
%                       disper(10) = ql phase factor (multiples of pi)
% 			disper(11) = Lorentzian width in energy.
% 			disper(12) = Temperature (K)
%
% Output variables:     s = intensity
%
%

%----- Calculate Dispersion

wq=sqrt(disper(4)^2+(disper(1)^2-disper(4)^2)*sin(disper(8) *pi*(qh-disper(5))).^2 +...
                    (disper(2)^2-disper(4)^2)*sin(disper(9) *pi*(qk-disper(6))).^2 +...
                    (disper(3)^2-disper(4)^2)*sin(disper(10)*pi*(ql-disper(7))).^2);

%----- Temperature factor

n_w=1./(exp(w*11.609/disper(12))-1);

%----- Calculate Intensity 

s=disper(11)^2*w.*(n_w+1).*...
    (1./((w-wq).^2+disper(11)^2)+1./((w+wq).^2+disper(11)^2));                                     
