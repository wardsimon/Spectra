function sout = exp(varargin)
%
% function r = log(s1..sn)
%
% @SPEC1D/exp function to give the natural exponential of each spectrum s1...sn.
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

s1 = [varargin{:}];

for i = 1:length(s1)
    r = s1(i);
    r.x = s1(i).x;
    r.y = exp(s1(i).y);
    r.e = s1(i).e.*r.y;
    if ~isempty(s1(i).yfit)
        r.yfit = exp(s1(i).yfit);
    end
    sout(i)=spec1d(r);
end