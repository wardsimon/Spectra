function conv=K2mev(dat)
% 1mev  =   1e-3 ev
% 1ev   =   1.6e-19 J
% dat(in J)/Kb= T(K)

dat=dat/1e-3;
dat=dat/1.60217653e-19;
conv=dat*1.3806504e-23;