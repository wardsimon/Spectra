function [sout,fitdata]=fits_win(s1,func,pin,notfixed,barriers,fcp,win)
% CK modified may 2008 from fits.m


% [sout,fitdata]=fits(s1,func,pin,notfixed,fcp,win)
%
% @SPEC1D/FIT Fits data in spec1d object s1 to MFIT function 
% specified in func. Also works for an array of spec1d objects.
%
% Notes:
% 1. If pin='auto' (or is ommitted) only s1 and func are specified, 
%    and func is a simple peak shape,
%    then an attempt is made to guess the start parameters.
% 2. Retrun parameters are:
%       r - original spectrum with fitted curve 
%       fitdata - structure with fields:                
%          fitdata.pvals - parameter values        
%          fitdata.evals - error values
%          fitdata.pnames- parameter names
%          fitdata.chisq - Chi^2 of fit
%
% Version 2.0, February 2001
% Des McMorrow and Henrik Ronnow

%--- Set defaults

if nargin==2, pin='auto'; notfixed=ones(1,4); end
if nargin==3, notfixed=ones(size(pin)); end
if ~exist('notfixed'), notfixed=ones(size(pin)); end
if nargin<=5, fcp=[0.0001 50 0.0001]; end
if isempty(fcp), fcp=[0.0001 50 0.0001]; end

if isempty(pin), pin='auto'; end 
if isempty(notfixed), notfixed=ones(1,4); end

%----- Loop over the number of spec1d objects

sout=spec1d;
fitdata=[];
pinin=pin;
for il=1:length(s1)

%----- Check for auto guess for simple peaks

   peaklist=strvcat('gauss','gauss2','lorz','lorz2','gauss_area');
   if ischar(pinin) %| nargin==2
   if ~isempty(strmatch(func,peaklist))
      stats=peakm(s1(il),1);      
      pin=[stats(1) stats(2) stats(3) stats(5)];   
   end   
   end   

   x=s1(il).x; y=s1(il).y; e=s1(il).e;    
    
%----- Check for window to remove points

   if nargin==6
     in=[];
     for n=1:size(win,1)
       in=[in x>win(n,1) & x<win(n,2)];  
     end
     x=x(in);
     y=y(in); 
     e=e(in);    
   end

%----- Remove zeros from e
   ezeros=find(e==0);
   x(ezeros)=[]; y(ezeros)=[]; e(ezeros)=[];

%----- Fit data
   [p, sig]=speclsqr_win(x,y,e,pin,notfixed,func,barriers,fcp);

%----- Fit function values
   yfit=feval(func,x,p);

%----- Calculate Chi^2
   v=length(y)-length(find(notfixed));
   ChiSq = sum( ((y-yfit)./e).^2 )/v;
   yfit=feval(func,x,p);

   ChiSq
   p'
%----- Get names of fit variables
   [dummy,dummy,pnames]=feval(func,x,p,1);

%----- set return
   sloop.x=x; sloop.y=y; sloop.e=e;
   sloop.x_label=s1(il).x_label; sloop.y_label=s1(il).y_label; sloop.datafile=s1(il).datafile;
   sloop.yfit=yfit;

   sout(il)=spec1d(sloop);
   
   fitdata(il).pvals=p;
   fitdata(il).evals=sig;
   fitdata(il).pnames=pnames;
   fitdata(il).chisq=ChiSq;

end
