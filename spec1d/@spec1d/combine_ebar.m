function r=combine(varargin)
%
% function r=combine([method,toll],s1,s2,....sn)
%
% @SPEC1D/COMBINE function to combine two or more spectra. 
%
% If the x values of two points differ by less
% than tolerance toll, then the points are combined.
%
% Depending on method, points are combined as
% 'mean'		: Simple means for x and y, errors are averaged in quadrature.
% 'counts'	: Restablishes normalisation and original counts assuming 
%				  square-root statistics, is correct for normalised counts
% 'weight'	: Weights to inverse error. For more general data.
% Default is 'counts'
%
% s1,s2,... can be single spectra or arrays of spectra.
%
% Example: 
% Combine s1,s2 and s3 when x values differ by less than 0.01.
% >r=combine('mean',0.01,s1,s2,s3)
%
% DFM 1.4.98, HMR, NBC, BHL 20.11.2000
x=[];
y=[];
e=[];
m=[];

% If no method specified, use 'counts' and shift input accordingly.
nextargin=1;
if isa(varargin{1},'char')
   method=varargin{1};
   nextargin=nextargin+1;
else
   method='counts';
end
if ~(strcmp(method,'mean') | strcmp(method,'counts') | strcmp(method,'weight'))
   warning('Unrecognized method, using ''counts''')
   method='counts';
end
if isa(varargin{nextargin},'double')
   toll=varargin{nextargin};
   nextargin=nextargin+1;
   notoll=0;
else
   notoll=1;
end

if notoll
   toll=min(diff(varargin{nextargin}(1).x))/10;
   if isempty(toll)
      toll=eps;
   end
end

% Gather all data points in one array.
for i=nextargin:length(varargin)
   if ~isa(varargin{i},'spec1d')
      disp('Append error: all objects must be spectra')      
      r=[];
      return
   end
   for n=1:length(varargin{i})
     if notoll & length(varargin{i}(n).x)>1
       toll=min(toll,min(diff(varargin{i}(n).x)/10));
     end
     x=[x; varargin{i}(n).x];
     y=[y; varargin{i}(n).y];
     e=[e; varargin{i}(n).e];
     % Reconstruct the montor count.
     % For zero error, this is impossible, so we take average monitor value.
     % This should really be replaced by a version keeping track of the monitor.
     mnew=zeros(size(varargin{i}(n).x));
     ne0=(varargin{i}(n).e~=0);
     mnew(ne0)=varargin{i}(n).y(ne0)./varargin{i}(n).e(ne0).^2;
     if any(ne0)
        mnew(~ne0)=mean(mnew(ne0));
     else
        mnew(~ne0)=mean(m);
     end
     m=[m; mnew];
   end
end
[x,n]=sort(x);
y=y(n);
e=e(n);
m=m(n);

if strcmp(method,'counts') & any(y<0)
   warning('Negative y-values. Only use combine_counts for counting data.')
end

% Take into account if error is zero
if strcmp(method,'weight') 
   x(e==0)=[];
   y(e==0)=[];
   e(e==0)=[];
end

xres=[];  % will contain final combined values
yres=[];
eres=[];
xcombi=[];   % will store data to be combined
ycombi=[];
ecombi=[];
mcombi=[];
for i=1:length(x);
   if isempty(xcombi)	
      xcombi=[x(i)];
      ycombi=[y(i)];
      ecombi=[e(i)];
      mcombi=[m(i)];
   else
      if 1/sqrt(sum(ycombi.*mcombi))<toll % keep combining untill relative error small enough
          if strcmp(method,'mean') % Simple mean
            xres=[xres;sum(xcombi)/length(xcombi)];
            yres=[yres;sum(ycombi)/length(ycombi)];
            eres=[eres;sqrt(sum(ecombi.*ecombi))/length(ecombi)];
         elseif strcmp(method,'counts') % Normalised counts
            xres=[xres;     sum(xcombi.*mcombi) /max(sum(mcombi),eps)];
            yres=[yres;     sum(ycombi.*mcombi) /max(sum(mcombi),eps)];
            eres=[eres;sqrt(sum(ycombi.*mcombi))/max(sum(mcombi),eps)];
         elseif strcmp(method,'weight') % Weighted by inverse error
            xres=[xres; sum(xcombi./ecombi)/sum(1./ecombi)];
            yres=[yres; sum(ycombi./ecombi)/sum(1./ecombi)];
            eres=[eres;sqrt(length(ecombi))/sum(1./ecombi)];
         end
         xcombi=[x(i)];
         ycombi=[y(i)];
         ecombi=[e(i)];
         mcombi=[m(i)];
      else
         xcombi=[xcombi;x(i)];
         ycombi=[ycombi;y(i)];       
         ecombi=[ecombi;e(i)];		
         mcombi=[mcombi;m(i)];		
      end
   end
end;

if ~isempty(xcombi)
   if strcmp(method,'mean') % Simple mean
      xres=[xres;sum(xcombi)/length(xcombi)];
      yres=[yres;sum(ycombi)/length(ycombi)];
      eres=[eres;sqrt(sum(ecombi.*ecombi))/length(ecombi)];
   elseif strcmp(method,'counts') % Normalised counts
      xres=[xres;     sum(xcombi.*mcombi) /max(sum(mcombi),eps)];
      yres=[yres;     sum(ycombi.*mcombi) /max(sum(mcombi),eps)];
      eres=[eres;sqrt(sum(ycombi.*mcombi))/max(sum(mcombi),eps)];
   elseif strcmp(method,'weight') % Weighted by inverse error
      xres=[xres; sum(xcombi./ecombi)/sum(1./ecombi)];
      yres=[yres; sum(ycombi./ecombi)/sum(1./ecombi)];
      eres=[eres;sqrt(length(ecombi))/sum(1./ecombi)];
   end
end;

%--- Create return spectrum
r.x=xres;
r.y=yres;
r.e=eres;
r.x_label=varargin{nextargin}(1).x_label;
r.y_label=varargin{nextargin}(1).y_label;
r.datafile=[];
r.yfit=[];

r=spec1d(r);
