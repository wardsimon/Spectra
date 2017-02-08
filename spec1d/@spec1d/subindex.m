function s_out = subindex(obj,ind)
for i = 1:length(obj)
    r = obj(i);
    r.x = r.x(ind);
    r.y = r.y(ind);
    r.e = r.e(ind);
    if ~isempty(r.yfit)
        r.yfit = r.yfit(ind);
    end
    s_out(i) = spec1d(r);
end
end