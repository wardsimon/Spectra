function [ varargout ] = feval(fit,varargin)
%
% function r = feval(fit,s1...sn)
%
% @SPEC1D/feval function to evaluate the results of a fit for one or more
% spectra.
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

p = inputParser;
p.CaseSensitive = false;
p.addRequired('fit',@(x) isstruct(x) & isfield(x,'function'));
p.addRequired('s_in',@(x) isa(x,'spec1d'));
p.addOptional('add_s',[],@(x) isa(x,'spec1d'))
p.parse(fit,varargin{:});

if ~isempty(p.Results.add_s(:))
    s = [p.Results.s_in(:) p.Results.add_s(:)];
else
    s = p.Results.s_in(:);
end


fit = p.Results.fit;

for i = 1:length(s)
    varargout{1}(i) = s(i);
    varargout{1}(i).yfit = feval(fit.function,s(i).x,fit.pvals);
end
