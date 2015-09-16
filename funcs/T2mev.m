function varargout=T2mev(varargin)
% varargout=T2mev(varargin)
% Converts Tesla to meV. The last argument is always g and we have assumed S = 1!
% -- Tesla * (S*g*ub) =  meV --
% Simon Ward Sept 2012 simon.ward@psi.ch

uB=5.7883817555E-2;

if nargin==2
    varargout{1}=uB*varargin{2}*varargin{1};
elseif nargin > 2
    for i=1:(nargin-1)
        varargout{i}=varargin{i}*(uB*varargin{nargin});
    end
end