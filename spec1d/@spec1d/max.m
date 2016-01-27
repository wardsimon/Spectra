function varargout = max(varargin)
%
% function [m ind xval] = max(s1..sn)
%
% @SPEC1D/max function to give the maximum of each spectrum s1...sn
% and optionally the index in the spectrum. Also the x-value of the corrsponding
% maximum can be returned.
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

s = [varargin{:}];

for i=1:length(s)
    [y_max, ind]=max(s(i).y);
    if nargout==1
        varargout{1}(i)=y_max;
    elseif nargout<=2
        varargout{1}(i)=y_max;
        varargout{2}(i)=ind;
    elseif nargout<=3
        varargout{1}(i)=y_max;
        varargout{2}(i)=ind;
        varargout{3}(i)=s(i).x(ind);
    end
end


