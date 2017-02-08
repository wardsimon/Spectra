function r=transform(transform_function,s1,varargin)
%
% function r=transform(transform_function,s1,params)
%
% @SPEC1D/TRANSFORM function to transform spectrum  s1 according to
% the file transform_function
%
% Note: varargin is a cell array that contains optional parameters
%
% Valid transforms: 1. normalise
%                   2. dydx
%
% Example
%
% Normalise spectrum s1
% >r=transform('normalise',s1,1);
% Differentiate spectrum s1
% >r=transform('dydx',s1);
%
% DFM 1.4.98
%
r=[];

if ~isa(s1,'spec1d')

   disp('Transform error: input spectrum must be a valid spec1d object')
   return

end

params=varargin;
[xt,yt,et,yfitt]=feval(transform_function,s1,params);

if isempty(xt); return; end
	
r.x=xt;
r.y=yt;
r.e=et;
r.x_label=s1.x_label;
r.y_label=s1.y_label;
r.datafile=[];
r.yfit=yfitt;
r=spec1d(r);
