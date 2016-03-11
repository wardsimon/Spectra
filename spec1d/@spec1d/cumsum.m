function s_out = cumsum(varargin)
%
% function r = cumsum(s1..sn)
%
% @SPEC1D/cumsum function to give the cumulative summation of points in each spectrum s1...sn.
%
% Simon Ward 27/01/2016 - simon.ward@psi.ch
%

s_ind = cellfun(@(x) isa(x,'spec1d'),varargin);
s = varargin(s_ind);
varargin = varargin(~s_ind);
s = [s{:}];

for i = 1:length(s);
    r = s(i);
    if isempty(varargin)
        r.y = cumsum(s(i).y);
    else
        r.y = cumsum(s(i).y,varargin{:});
    end
    r.e = sqrt(cumsum(s(i).e.^2));
    
    if ~isempty(s(i).yfit)
        if isempty(varargin)
            r.yfit = cumsum(s(i).yfit);
        else
            r.yfit = cumsum(s(i).yfit,varargin{:});
        end
    end
    
    s_out(i) = feval(class(r),r);
end