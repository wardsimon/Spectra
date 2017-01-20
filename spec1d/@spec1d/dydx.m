function r=dydx(s1)
%
% function r=dydx(s1)
%
% @SPEC1D/DYDX function to differentiate spectrum s1.
%
% DFM 1.4.98 
%
yfitd=[];

x=s1.x; x=x(:);
y=s1.y; y=y(:);
e=s1.e; e=e(:);
yfit=s1.yfit; yfit=yfit(:);

yd=diff(y)./diff(x);
ys=[NaN;yd];
ys(end)=[];
yd=yd+ys;
yd(1)=[];

if ~isempty(yfit)
  yfitd=diff(yfit)./diff(x);
  ys=[NaN;yfitd];
  ys(end)=[];
  yfitd=yfitd+ys;
  yfitd(1)=[];
end  
  
yd=0.5*yd;
xd=x(2:length(x)-1);
ed=e(2:length(e)-1);
yfitd=0.5*yfitd;

r.x=xd;
r.y=yd;
r.e=ed;
r.x_label=s1.x_label;
r.y_label=s1.y_label;
r.datafile=s1.datafile;
r.yfit=yfitd;

r=spec1d(r);