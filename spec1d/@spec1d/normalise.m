function sn=normalise(s1,normval)
%
% @SPEC1D/NORMALISE function to normalise spectrum to
% value specified in params.
%
% DFM 1.4.98
% Modified by JOP 20.6.10

if nargin < 2
    normval=1;
end

for n=1:length(s1)
    x=s1(n).x;
    y=s1(n).y;
    e=s1(n).e;
    yfit=s1(n).yfit;

    xt=x;
    ymax=max(y);
    yt=y/ymax*normval;
    et=e/ymax*normval;
    yfitt=yfit/ymax*normval;

    r.x=xt;
    r.y=yt;
    r.e=et;
    r.x_label=s1(n).x_label;
    r.y_label=s1(n).y_label;
    r.datafile=s1(n).datafile;
    r.yfit=yfitt;

    sn(n)=spec1d(r);
end