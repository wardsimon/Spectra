function r1=log(s1)
%
% function r=abs(s1)
%
% @SPEC1D/DYDX function to give the absolute value of spectrum s1.
%
% JOP 4.11.09
%
yfitabs=[];

for n=1:length(s1)
    x=s1(n).x; x=x(:);
    y=s1(n).y; y=y(:);  yy=y(:);
    e=s1(n).e; e=e(:);
    yfit=s1(n).yfit; yfit=yfit(:);
    s = warning('query', 'all');
    warning off all
    
    yabs=log(y);

    if ~isempty(yfit)
        yfitabs=log(yfit);
    end
    yabs(yy<0)=0;
    
    r.x=x;
    r.y=yabs;
    r.e=log(e)./e;
    r.e(yy<0)=0;
    warning(s);
    r.x_label=s1(n).x_label;
    r.y_label=s1(n).y_label;
    r.datafile=s1(n).datafile;
    r.yfit=yfitabs;

    r1(n)=spec1d(r);
end