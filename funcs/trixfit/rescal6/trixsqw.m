function [s]=trixsqw(disper,qh,qk,ql,w)
%
% MATLAB  function to define S(Q,w) for TRIX
%
% Des Mcmorrow 7.11.95
% 
% Notes: (1) The lineshape is assumed to be Lorentzian
%        (2) S(Q,w) is assumed to follow a 3D dispersion
%        (3) I have not added the kf/ki term yet!
%
% Called by: TRIX
% Calls  to: 
%
% Units: At the moment will only work with meV
%
% Input variables:      disp = dispersion parameters
%			disper(1) = qh zone-boundary energy
%			disper(2) = qk      "          "
% 			disper(3) = ql      "          "
% 			disper(4) = zone-centre energy gap
% 			disper(5) = qh zone-centre
% 			disper(6) = qk zone-centre
% 			disper(7) = ql zone-centre
% 			disper(8) = Lorentzian width in energy.
% 			disper(9) = Temperature (K)
%
% Output variables:
%

%----- Calculate Dispersion

wq=sqrt(disper(4)^2+(disper(1)^2-disper(4)^2)*sin(pi*(qh-disper(5))).^2 +...
                    (disper(2)^2-disper(4)^2)*sin(pi*(qk-disper(6))).^2 +...
                    (disper(3)^2-disper(4)^2)*sin(pi*(ql-disper(7))).^2);

%----- Temperature factor

n_w=1./(exp(w*11.609/disper(9))-1);

%----- Calculate Intensity 

%----- Damped harmonic oscillator
%
% s=disper(8)*w./q.*(n_w+1)./((w.^2-wq.^2).^2+(disper(8)*w).^2);

%----- Sum of Lorentzians
                                 
s=disper(8)^2*w.*(n_w+1).*(1./((w-wq).^2+disper(8)^2)+1./((w+wq).^2+disper(8)^2));  

%----- 1/q factor for antiferromagnets
%
% q=sqrt((qh-disper(5)).^2+(qk-disper(6)).^2+(ql-disper(7)).^2);
% s=s./q;

