function [varargout]=mapplot(s,varargin)
    %%
    % FUNCTION  : [xi,yi,zi]=mapplot(varargin)
    % VARARGIN  : Required:
    %                 s       = spec1d objects, in the form s(1:2) or [s1 s2]
    %             Optional:
    %                 spacing = Vector of length s. 1:length(s) will be used otherwise
    %                 smooth  = Gaussian FWHM to smooth data. length(smooth)=1;
    %                 x_bin   = Vector for x binning. Note length(x_bin) ~= 1 or length(s)
    %                 y_bin   = Vector for y binning. Note length(y_bin) ~= 1 or length(s)
    %                               !! y_bin must immediately follow x_bin !!
    %               pre_fun = Functions to be applied to z data. i.e log, real etc
    % VARARGOUT : No outputs plots, otherwise [x_data, y_data, z_data] 
    % NOTES     : If you want to bin in x but not y, you need to supply y_bin=[].
    %             i.e. mapplot(... x_bin,[],....)
    %             WE DO NOT CHECK THE VALIDITY OF pre_fun!
    % AUTHOR    : Simon Ward
    % DATE      : Version 1, 23/04/2013
    
    extvar=1:length(s);
    xi=[];
    yi=[];
    smoothe=[];
    scale = {};
    j=1;
    f_x=0;
    for i=1:length(varargin)
        temp=varargin{i};
        if isnumeric(temp)
            if length(temp)==1
                smoothe=temp;
            else
                if length(temp)==(length(s))
                    extvar=temp;
                else
                    if ~f_x
                        xi=temp;
                        yi=varargin{i+1};
                        f_x=1;
                    else
                        yi=temp;
                    end
                end
            end
        elseif ischar(temp)
            scale{j}=temp;
            j=j+1;
        elseif isempty(temp)
            % f_x is 1 if xi = [];
            f_x=1;
        end
    end
    
    x=[];
    y=[];
    z=[];
    for n=1:length(s)
        if ~isempty(smoothe)
            s(n)=smooth(s(n),smoothe);
        end
        y=[y;s(n).('x')];
        z=[z;s(n).('y')];
        if size(extvar)==[length(s(n).('x')) length(s)]
            x=[x;extvar(:,n)];
        else
            x=[x;extvar(n)*ones(size(s(n).('x')))];
        end
    end
    
    if ~isempty(scale)
        for i=1:length(scale)
            temp=scale{i};
            f=str2func(temp);
            z=feval(f,z);
        end
    end
    
    if isempty(xi)
        xi = linspace(min(x),max(x),200);
    end
    if isempty(yi)
        yi = linspace(min(y),max(y),200);
    end
    xi = xi(:)';
    yi = yi(:)';
    
    
    if exist('scatteredInterpolant','file')==2
        F = scatteredInterpolant(x(:),y(:),z(:),'natural');
        [xi, yi] = meshgrid(xi,yi);
        zi=F(xi,yi);
    else
        zi = griddata(x,y,z,xi,yi','cubic');
        [xi, yi] = meshgrid(xi,yi);
    end
    
    if nargout == 0
        figure
        p = pcolor(xi,yi,zi);
        set(p,'LineStyle','None')
        shading flat
        line(x,y,'linestyle','none','marker','.','color','k')
    elseif nargout == 3
        varargout{1} = xi;
        varargout{2} = yi;
        varargout{3} = zi;
    end
end
