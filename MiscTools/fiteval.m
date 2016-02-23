function [ varargout ] = feval(fit,varargin)
%
% function r = feval(fit,s1...sn)
%
% @SPEC1D/feval function to evaluate the results of a fit for one or more
% spectra.
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

% Deal with the case of someone not knowing seval.
s_ind = cellfun(@(x) isa(x,'spec1d'),varargin);
s = varargin(s_ind);
varargin(s_ind) = [];

if ~isempty(s)
    varargout = cell(size(s));
   [varargout{:}] = seval(fit,s); 
   return
end

p = inputParser;
p.CaseSensitive = false;
p.addRequired('fit',@(x) isstruct(x) & isfield(x,'function'));
p.addRequired('x',@isnumeric);
p.parse(fit,varargin{:});

x = [p.Results.x];
fit = p.Results.fit;

for i = 1:length(fit)
    if length(fit) == nargout
        varargout{i}(:,i) = feval(fit.function,x,fit.pvals);
    else
    varargout{1}(:,i) = feval(fit.function,x,fit.pvals);
    end
end
