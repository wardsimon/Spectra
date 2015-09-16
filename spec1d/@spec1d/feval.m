function [ varargout ] = feval( s,fit )
    %EVAL Summary of this function goes here
    %   Detailed explanation goes here
    for i = 1:length(s)
        varargout{1}(i) = s(i);
        varargout{1}(i).yfit = feval(fit.function,s(i).x,fit.pvals);
    end
end

