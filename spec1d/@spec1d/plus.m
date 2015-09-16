function rvec=plus(s1vec,s2)
    %
    % function r=plus(s1,s2)
    %
    % @SPEC1D/PLUS plus for 1D spectra
    %
    % Valid expressions:
    %
    % 1. r=s1+s2     , where s1 and s2 are spectra
    % 2. r=s1+a      , where a is to be added to the y axis
    % 3. r=s1+[xc,yc], where xc (yc) is added the x (y) axis
    %
    % If spectra are not the same length, s2 is interpolated to
    % the x values of s1.
    %
    % Plus should not be used to combine counts taken with different
    % monitors - use combine(s1,s2)
    %
    % DFM 1.4.98 HMR 20.11.2000
        
    for nvec=1:length(s1vec)
        s1=s1vec(nvec);
        
        r=[];
        xs=0;
        yfit=[];
        yfitch=0;
        
        if isa(s1,'spec1d') & isa(s2,'spec1d')
            if length(s1.x)~=length(s2.x)
                disp('Spectra are not the same length')
                disp('Uses interpolation for second spectrum')
                s2=interpolate(s2,s1);
            elseif (s1.x-s2.x) > mean(diff(s1.x))/10
                disp('X-values not the same')
                disp('Uses interpolation for second spectrum')
                s2=interpolate(s2,s1);
            end
            x1=s1.x; y1=s1.y; e1=s1.e;
            x2=s2.x; y2=s2.y; e2=s2.e;
            x_label=s1.x_label; y_label=s1.y_label;
            if any(abs(log(abs((y1./e1.^2)./(y2./e2.^2))))>0.0953)
                warning('Plus should not be used to add counts with different monitors - use combine(s1,s2)');
            end
        elseif ~isa(s2,'spec1d')
            x1=s1.x; y1=s1.y; e1=s1.e; yfit=s1.yfit;
            if length(s2)==2
                x2=x1; xs=s2(1); y2=s2(2); yfitch=s2(2);
            elseif length(s2)==1;
                x2=x1; y2=s2(1); yfitch=s2(1);
            elseif length(s2)==length(y1)
                x2=x1;
                y2=s2;
                yfitch=s2;
            end
            e2=0;
            x_label=s1.x_label; y_label=s1.y_label;
        elseif ~isa(s1,'spec1d')
            x2=s2.x; y2=s2.y; e2=s2.e; yfit=s2.yfit;
            if length(s1)==2
                x1=x2; xs=s1(1); y1=s1(2); yfitch=s1(2);
            elseif length(s1)==1;
                x1=x2; y1=s1(1); yfitch=s1(1);
            elseif length(s1)==length(y2)
                x1=x2; y1=s1;  yfitch=s1;
            end
            e1=0;
            x_label=s2.x_label; y_label=s2.y_label;
        end
        
        r.x=x1+xs;
        r.y=y1+y2;
        r.e=sqrt(e1.^2+e2.^2);
        r.x_label=x_label;
        r.y_label=y_label;
        r.datafile=[];
        if length(yfit)==length(yfitch)
            r.yfit=yfit+yfitch;
        end
        r=spec1d(r);
        
        rvec(nvec)=r;
        
    end
    
    return