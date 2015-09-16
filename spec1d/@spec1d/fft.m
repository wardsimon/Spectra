function sout=fft(varargin)
k=1;
for i=1:length(varargin)
    for j=1:length(varargin{i})
        if isa(varargin{i}(j),'spec1d')
            s1(k)=varargin{i}(j);
            k=k+1;
        end
    end
end

for i=1:(k-1)
    x=s1(i).x;
    y=s1(i).y;
    
    L=length(x);
    NFFT = 2^nextpow2(L);
    Y = fft(y,NFFT)/L;
    dx=diff(x)';
    f = L*(dx(1:NFFT/2+1)/2).*linspace(0,1,NFFT/2+1);
    Y=Y(1:NFFT/2+1);
    E=fft(s1(i).e,NFFT)/L;
    E=E(1:NFFT/2+1);
    
    r.x=f;
    r.y=Y;
    r.e=E;
    if ~isempty(s1(i).yfit)
        yf=fft(s1(i).yfit,NFFT)/L;
        r.yfit=yf(1:NFFT/2+1);
    end
    r.x_label=s1(i).x_label;
    r.y_label=s1(i).y_label;
    r.datafile=s1(i).datafile;
    sout(i)=spec1d(r);
end