function [varargout]=differentiate(varargin)

j=1;
for n=1:length(varargin)
     for i=1:length(varargin{n})
        x=varargin{n}(i).x;
        y=varargin{n}(i).y;
        e=varargin{n}(i).e;
        xn=x;   yn=y;   en=e;
        for k=1:nargout
            xn=xn(1:end-1)+diff(xn)/2;
            yn=diff(yn);
            en=en(1:end-1);
            varargout{k}(j)=spec1d(xn,yn,en);
            j=j+1;
        end
     end
%    f=fittype('pchipinterp');
%    Linfit=fit(x,y,f);
%    [d1(n,:) d2(n,:)]=differentiate(Linfit,x);
%    svec0(n)=spec1d(x,d1(n,:),e);
%    svec1(n)=spec1d(x,d2(n,:),e);
end