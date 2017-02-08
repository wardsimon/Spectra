function sout=sign(varargin)

s1 = cellfun(@(x)isa(x,)))

for i=1:(k-1)
    r = s1(i);
    r.x = s1(i).x;
    r.y = sign(s1(i).y);
    r.e = zeros(size(s1(i).x));
    if ~isempty(s1(i).yfit)
        r.yfit = sign(s1(i).yfit);
    end
    r.x_label=s1(i).x_label;
    r.y_label=s1(i).y_label;
    r.datafile=s1(i).datafile;
    sout(i) = feval(class(r),r);
end