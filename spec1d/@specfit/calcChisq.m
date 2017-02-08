function calcChisq( obj,s )
%CALCCHISQ Summary of this function goes here
%   Detailed explanation goes here
if ~exist('s','var')
    try
        s = obj.findspec1d();
    catch
        
    end
end
        % Goodness of fit
        v = length(s.y)-sum(logical(obj.notfixed));
        obj.chisq = sum(((s.y-s.yfit)./s.e).^2 )/v;
%         ret = obj.chisq;

end

