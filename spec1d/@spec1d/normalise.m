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

for n = 1:length(s1)
    x = s1(n).x;
    y = s1(n).y;
    e = s1(n).e;
    if ~isempty(s1(n).yfit)
        yfit = s1(n).yfit;
    end
    
    ymax = max(y);
    yt = y/ymax*normval;
    et = e/ymax*normval;
    if ~isempty(s1(n).yfit)
        yfitt = yfit/ymax*normval;
    else
        yfitt = [];
    end
    
    r = s1(n);
    r.y = yt;
    r.e = et;
    r.yfit = yfitt;

    sn(n) = feval(class(r),r);
end