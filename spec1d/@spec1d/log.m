function sout = log(varargin)
%
% function r = log(s1..sn)
%
% @SPEC1D/Log function to give the natural log of each spectrum s1...sn.
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

s1 = [varargin{:}];

for n = 1:length(s1)
    
    x = s1(n).x(:);
    y = s1(n).y(:);
    e = s1(n).e(:);
    yfit = s1(n).yfit(:);
    
    r = s1(n);
    r.x = x;
    r.y = log(y);
    r.e = e./y;
    
    if ~isempty(yfit)
        r.yfit = log(yfit);
    end
    sout(n) = feval(class(r),r);
end