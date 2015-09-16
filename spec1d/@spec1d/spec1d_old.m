function s=spec1d(a,b,c)
%
% function s=spec1d(a)
%
% MATLAB function to create a spec1d object s
%
% DFM 1.4.98 
%
s.x=[];
s.y=[];
s.e=[];
s.x_label='';
s.y_label='';
s.datafile='';
s.yfit=[];
if nargin==0
   s=class(s,'spec1d');
elseif nargin>2
%    [a,n]=sort(a);
   s.x=a(:);
   if length(b)==1
      b=b+0*s.x;
   end
%    b=b(n);
   s.y=b(:);
   if length(c)==1
      c=c+0*s.x;
   end
%    c=c(n);
   s.e=c(:);
   s=class(s,'spec1d');
elseif isa(a,'spec1d')
   s=a;
elseif isa(a,'psd_obj')
   a=struct(a);
   nspec=size(a.WinSum,2);
   for n=1:nspec
     s(n).x=a.X(:);
     s(n).y=a.WinSum(:,n);
     s(n).e=a.WinSumErr(:,n);
     s(n).x_label=a.XLab;
     s(n).y_label='PSD data';
     s(n).datafile=a.FileName;
     s(n).yfit=[];
   end
   s=class(s,'spec1d');
elseif isa(a,'struct')
   s.x=a.x(:);
   s.y=a.y(:);
   s.e=a.e(:);
   if isfield(a,'x_label')
      s.x_label=a.x_label;
   else
      s.x_label='';
   end
   if isfield(a,'y_label')
      s.y_label=a.y_label;
   else
      s.y_label='';
   end
   if isfield(a,'datafile')
      s.datafile=a.datafile;
   else
      s.datafile='';
   end
   if isfield(a,'yfit')
      s.yfit=a.yfit;
   else
      s.yfit=[];
   end
   s=class(s,'spec1d');
elseif isa(a,'double')
   if size(a,2)==3
     [aa,n]=sort(a(:,1));
     s.x=a(:,1);
     s.y=a(:,2);
     s.e=a(:,3);
     s=class(s,'spec1d'); 
   else
      error('Cannot convert input to SPEC1D object')
   end
else
   error('Cannot convert input to SPEC1D object')
end

if any(diff([length(s.x) length(s.y) length(s.e)]))
    warning('x, y and e are not of the same length!\nLength x:\t%i\nLength y:\t%i\nLength e:\t%i',length(s.x),length(s.y),length(s.e))
end

% [s.x, ind ] = sort(s.x);
% s.y = s.y(ind);
% s.e = s.e(ind); 
