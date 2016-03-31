function s_out = loadobj(s_in)
%% This function fixes the conversion between the original spec1d and the post 2013 version.
s = size(s_in);

s_in = s_in(:);

for i = 1:length(s_in)
    s_in(i).x = s_in(i).x(:);
    [s_in(i).x, ind] = sort(s_in(i).x);
    s_in(i).y = s_in(i).y(:);
    s_in(i).y = s_in(i).y(ind);
    s_in(i).e = s_in(i).e(:);
    s_in(i).e = s_in(i).e(ind);
    if ~isempty(s_in(i).yfit)
        s_in(i).yfit = s_in(i).yfit(:);
        s_in(i).yfit = s_in(i).yfit(ind);
    end
    s_out(i) = spec1d(s_in(i));
end

s_out = reshape(s_out,s);
