function varargout = sum(varargin)
%
% function r = sum(s1...sn)
%
% @SPEC1D/sum function to give the sum of elements in each spectrum s1...sn.
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

s = [varargin{:}];

out_y = zeros(size(s));
out_e = out_y;

for i = 1:length(s)
    out_y(i) = sum(s(i).y);
    out_e(i) = sqrt(sum(s(i).e.^2));
end

varargout{1} = out_y;

if nargout == 2
    varargout{2} = out_e;
end