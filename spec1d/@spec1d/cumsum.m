function sout = cumsum(varargin)
%
% function r = cumsum(s1..sn)
%
% @SPEC1D/cumsum function to give the cumulative summation of points in each spectrum s1...sn.
%
% Simon Ward 27/01/2016 - simon.ward@psi.ch
%

s1 = [varargin{:}];

for i = 1:length(s1);
    r = s1(i);
    r.y = cumsum(s1(i).y);
    r.e = sqrt(sum(s(i).e.^2));
    
    if ~isempty(s1(i).yfit)
        r.yfit = cumsum(s1(i).yfit);
    end
    
    sout(i) = spec1d(r);
end