function r=rdivide(s1,s2)
%
% function r=rdivide(s1,s2)
%
% @SPEC1D/RDIVIDE Rdivide for 1D spectra
%
% Valid expressions: 
%
% 1. r=s1./s2      , where s1 and s2 are spectra
% 2. r=s1./yc      , where yc is a constant that divides the y-axis
% 3. r=s1./[xc,yc] , where xc (yc) divides the x (y) axis by a constant
%
% DFM 1.4.98
%
r=[];
xs=1;
yfit=[];
yfitch=0;

if isa(s1,'spec1d') & isa(s2,'spec1d')
  if length(s1.x)~=length(s2.x)
     disp('Spectra are not the same length')
     disp('Uses interpolation for second spectrum')
     s2=interpolate(s2,s1);
  end

  x1=s1.x; y1=s1.y; e1=s1.e;
  x2=s2.x; y2=s2.y; e2=s2.e;
  x_label=s1.x_label; y_label=s1.y_label;

elseif ~isa(s2,'spec1d')
    x1=s1.x; y1=s1.y; e1=s1.e; yfit=s1.yfit;
    if length(s2)==2
      x2=x1; xs=s2(1); y2=s2(2); yfitch=s2(2);
  elseif length(s2)==1;
     x2=x1; y2=s2(1); yfitch=s2(1);
  elseif length(s2)==length(x1)
      x2=x1; y2=s2;
  end
  e2=0;
  x_label=s1.x_label; y_label=s1.y_label;

elseif ~isa(s1,'spec1d')

  x2=s2.x; y2=s2.y; e2=s2.e; yfit=s2.yfit;
  if length(s1)==2
     x1=x2; xs=s1(1); y1=s1(2); yfitch=s1(2);
  elseif length(s1)==1;
     x1=x2; y1=s1(1); yfitch=s1(1);
  end
  e1=0;
  x_label=s2.x_label; y_label=s2.y_label;

end

if length(x1)~=length(x2)
   disp('Spectra are not the same length')
   return
end

r.x=x1./xs;
r.y=y1./y2;
r.e=sqrt((y1.*e2).^2+(y2.*e1).^2)./y2.^2;
r.x_label=x_label;
r.y_label=y_label;
r.datafile=[];
r.yfit=[];
if yfitch~=0 r.yfit=yfit/yfitch; end

r=spec1d(r);






