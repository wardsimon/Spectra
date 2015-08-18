function varargout=mev2K(varargin)
% varargout=mev2K(varargin)
% Converts meV to Kelvin.
% Simon Ward Sept 2012 simon.ward@psi.ch

for i=1:nargin
    % data in mJ/kB
    varargout{i}=varargin{i}*1e-3*1.60217653e-19/1.3806504e-23;
end