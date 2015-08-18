function [varargout]=mapplots(s,split,extvar,varargin)
    %%
    % FUNCTION  : [xi,yi,zi]=mapplots(s,split,extvar,varargin)
    % VARARGIN  : Required:
    %                 s       = spec1d objects, in the form s(1:2) or [s1 s2]
    %                 split   = where to split the object (remember if [s1 s2] it'll be length(s1)+1 !)
    %                 spacing = cell containing Vectors which must be of length s1 and s2 i.e. {1:length(s1);1:length(s2)}
    %             Optional:
    %                 smooth  = Gaussian FWHM to smooth data. length(smooth)=1;
    %                 x_bin   = Vector for x binning. Note length(x_bin) ~= 1 or length(s)
    %                 y_bin   = Vector for y binning. Note length(y_bin) ~= 1 or length(s)
    %                               !! y_bin must immediately follow x_bin !!
    %               pre_fun = Functions to be applied to z data. i.e log, real etc
    % VARARGOUT : No outputs plots, otherwise [x_data, y_data, z_data] 
    % NOTES     : If you want to bin in x but not y, you need to supply y_bin=[].
    %             i.e. mapplots(... x_bin,[],....)
    %             WE DO NOT CHECK THE VALIDITY OF pre_fun!
    % AUTHOR    : Simon Ward
    % DATE      : Version 1, 14/07/2015
    
    if split >= length(s)
        error('Split needs to be less than the length of spec1d objects')
    end
    
    s1 = s(1:(split)-1);
    s2 = s(split:end);
    
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
                if ~f_x
                    xi=temp;
                    yi=varargin{i+1};
                    f_x=1;
                else
                    yi=temp;
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
    
    for n=1:length(s1)
        if ~isempty(smoothe)
            s1(n)=smooth(s1(n),smoothe);
        end
    end
    for n=1:length(s2)
        if ~isempty(smoothe)
            s2(n)=smooth(s2(n),smoothe);
        end
    end
    
    [y1, z1, ~] = extract(s1);
    if iscell(y1)
        for n = 1:length(s1)
           x1{n} = extvar{1}(n)*ones(size(y1{n}));
        end
    else
        x1 = mat2cell(repmat(extvar{1}(:)',size(y1,1),1),length(y1),ones(size(s1)));
        y1 = mat2cell(y1,length(y1),ones(size(s1)));
        z1 = mat2cell(z1,length(z1),ones(size(s1)));
    end
    [x2, z2, ~] = extract(s2);
    if iscell(x2)
        for n = 1:length(s2)
           y2{n} = extvar{2}(n)*ones(size(x2{n}));
        end
    else
        y2 = mat2cell(repmat(extvar{2}(:)',size(x2,1),1),length(x2),ones(size(s2)));
        x2 = mat2cell(x2,length(x2),ones(size(s2)));
        z2 = mat2cell(z2,length(z2),ones(size(s2)));
    end
    x = [vertcat(x1{:});vertcat(x2{:})];
    y = [vertcat(y1{:});vertcat(y2{:})];
    z = [vertcat(z1{:});vertcat(z2{:})];
    
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
    
    
    if exist('TriScatteredInterp','file')==2
        F = TriScatteredInterp(x(:),y(:),z(:),'natural');
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
