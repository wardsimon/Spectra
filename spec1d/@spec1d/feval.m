function [ varargout ] = feval(fit,varargin)
%
% function r = feval(fit,s1...sn)
%
% @SPEC1D/feval function to evaluate the results of a fit for one or more
% spectra.
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%
s_ind = cellfun(@(x) isa(x,'spec1d'),varargin);
s = varargin(s_ind);
varargin(s_ind) = [];

p = inputParser;
p.CaseSensitive = false;
p.addRequired('fit',@(x) isstruct(x) & isfield(x,'function'));
p.addRequired('s_in',@iscell);
p.parse(fit,s,varargin{:});

s = [p.Results.s_in{:}];

fit = p.Results.fit;

for i = 1:length(s)
    varargout{1}(i) = s(i);
    varargout{1}(i).yfit = feval(fit.function,s(i).x,fit.pvals);
end
