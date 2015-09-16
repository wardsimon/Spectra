function rvec = times(varargin)
    s = {}; f = [];
    for i = 1:2
        if isa(varargin{i},'spec1d')
            s{end+1} = varargin{i};
        else
            f = varargin{i};
        end
    end
    
    if isempty(f)
        if length(s)~=2
            error('Somethings gone wrong.....')
        end
        % Multiply spec1d objects
        for i = 1:length(s{1})
            if length(s{1}(i).x)~=length(s{2}.x)
                warning('Objects are not the same length. Using interpolation for second object')
                s_2 = interpolate(s{2},s{1}(i));
            else
                s_2 = s{2};
            end
            rvec(i) = s{1}(i);
            rvec(i).y = rvec(i).y .* s_2.y;
            rvec(i).e = rvec(i).y.*sqrt((s{1}(i).e./s{1}(i).y).^2 + (s_2.e./s_2.y).^2);
            rvec(i).yfit = s{1}(i).yfit .* s_2.yfit;
        end
    else
        k = 1;
        for i = 1:length(s)
            for j = 1:length(s{i})
                rvec(k)   = s{i}(j);
                rvec(k).y = f.*rvec(k).y;
                rvec(k).e = f.*rvec(k).e;
                rvec(k).yfit = f.*rvec(k).yfit;
                k = k + 1; 
            end
        end
    end
    