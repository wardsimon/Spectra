function varargout = min(varargin)
%
% function r = min(s1..sn)
%
% @SPEC1D/min function to find the minimum value of each spectrum s1...sn.
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

s = [varargin{:}];

for i=1:length(s)
    [y_min, ind]=min(s(i).y);
    if nargout==1
        varargout{1}(i)=y_min;
    elseif nargout<=2
        varargout{1}(i)=y_min;
        varargout{2}(i)=ind;
    elseif nargout<=3
        varargout{1}(i)=y_min;
        varargout{2}(i)=ind;
        varargout{3}(i)=s(i).x(ind);
    end
end