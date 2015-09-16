function s_out = mpower(s,p)
    size_s = size(s);
    s = s(:);
    
    for i=1:length(s)
        s_out(i) = s(i);
        s_out(i).y = s(i).y.^p;
        s_out(i).e = p*(s(i).e./s(i).y).*s_out(i).y;
        if ~isempty(s(i).yfit)
            s_out(i).yfit = s_out(i).yfit.^p;
        end
    end
    
    s_out = reshape(s_out,size_s);