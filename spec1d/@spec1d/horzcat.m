function r=horzcat(varargin)
%
% function r=horzcat(s1,s1,...)
%
% SPEC1D/HORZCAT Creates an array of spectra
%
% Examples:
% 1. s3=[s1 s2], makes an array of spec1d objects s3 from s1 and s2
%
% HMR 30.11.2000;
m=1;
for i=1:length(varargin)

   if ~isa(varargin{i},'spec1d')
      disp('Append error: all objects must be spectra')      
      return
   end
   for n=1:length(varargin{i})
      r(m)=varargin{i}(n);
      m=m+1;
   end
end
