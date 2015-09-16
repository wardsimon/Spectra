function [p]=rc_savp2(pradat,cfgdat)
%
% MATLAB function rc_savp saves parameter  from RESCAL windows
%        either to a file or a
% Note: If option == 'pfile' , save RESCAL parameters to file *.par
%       If option == 'ifile' , save instrument parameters to file *.cfg
%       If option == 'sfile' , save simulation parameters to file *.sim
%       If option == 'spars' ,   unload simulation parametes to p
%       If option == 'respars' , unload all rescal parameters to p
%       If option == 'tpars'   , unload trix scan parameters to p
%
% DFM 9.5.95
%
cur_dir=cd;
pres = []; pinst = []; psim = []; ptrix = [];

pres=zeros([42,1]);

for i=1:length(pradat)
    pres(i)=pradat(i);
end
% pres=pradat(1:42);

pinst=zeros([27,1]);

pinst(1) =cfgdat(1);
pinst(2) =cfgdat(2);
pinst(3) =cfgdat(3);
if pinst(1)  == 1, pinst(3)=pinst(2); end              % Circular source

pinst(4) =cfgdat(4);
pinst(5) =cfgdat(5);
pinst(6) =cfgdat(6);
if pinst(4)  == 1, pinst(5)=0; pinst(6)=0; end         % No guide

pinst(7) =cfgdat(7);
pinst(8) =cfgdat(8);
pinst(9) =cfgdat(9);
pinst(10)=cfgdat(10);
if pinst(7)  == 1, pinst(10)=pinst(9); pinst(9)=pinst(8); end  % Cylindrical sample

pinst(11) =cfgdat(11);
pinst(12) =cfgdat(12);
pinst(13) =cfgdat(13);
if pinst(11) == 1, pinst(13)=pinst(12); end            % Circular detector

pinst(1)=pinst(1)-1;   pinst(4)=pinst(4)-1;   pinst(7)=pinst(7)-1;   pinst(11)=pinst(11)-1;


for i=14:length(cfgdat)
    pinst(i)=cfgdat(i);
end

p=[pres; pinst; psim; ptrix];
