function varargout=mev2T(varargin)
% varargout=mev2T(varargin)
% Converts meV to Tesla. The last argument is always g and we have assumed S = 1!
% -- mev/(S*g*ub) = Tesla --
% Simon Ward Sept 2012 simon.ward@psi.ch

uB=5.7883817555E-2;

if nargin==2
    varargout{1}=varargin{1}/(uB*varargin{2});
elseif nargin > 2
    for i=1:(nargin-1)
        varargout{i}=varargin{i}/(uB*varargin{nargin});
    end
end