function r=vertcat(varargin)
%
% function r=vertcat(s1,s2,...)
%
% SPEC1D/VERTCAT Appends spectra into single spec1d object
%
% Valid expressions: 
%
% 1. r=[s1; s2] , where s1 and s2 are spec1d objects
% 2. r=[s1; s2; s3], where s1, s2 and s3 are spec1d objects
%
% DFM 1.4.98;
%
r=varargin{1};

if length(varargin)>1
for i=2:length(varargin)

   if ~isa(varargin{i},'spec1d')
      disp('Append error: all objects must be spectra')      
      return
   end

   r.x=[r.x; varargin{i}.x];
   r.y=[r.y; varargin{i}.y];
   r.e=[r.e; varargin{i}.e];

end

  [r.x,srt]=sort(r.x);
  r.y=r.y(srt);
  r.e=r.e(srt);
  r.datafile=[];
  r.yfit=[];
  r=spec1d(r);
end


