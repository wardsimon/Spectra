function [y, name, pnames, pin]=fourier_series(x,p, flag)
% a0 + E ((aN*sin(2*pi*N*x) + bN*cos(2*pi*n*x))
    
order = (length(p)-2)/2;

if nargin==2
    y=p(1)*ones(size(x));
    
    w=p(2);
    
    for i=1:order
        off=2*i+1;
        aN=p(off);
        bN=p(off+1);
        y=y+(aN*sin((off-2)*x*w) + bN*cos((off-2)*x*w));
    end
else
    y=zeros(size(x));
    name=sprintf('%i Order Fourier Series',order);
    pnames={'a0','w'};
    for i=1:order
        off=2*i+1;
        pnames{off}=sprintf('a%i',i);
        pnames{off+1}=sprintf('b%i',i);
    end
    pin=[];
end
    
