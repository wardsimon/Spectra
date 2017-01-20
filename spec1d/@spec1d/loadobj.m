function s_out = loadobj(s_in)
%% This function fixes the conversion between the original spec1d and the post 2013 version.
sz = size(s_in);

s_in = s_in(:);

s = spec1d();
f = fieldnames(s);
s = repmat(s,size(s_in));

for i = 1:length(s_in)
<<<<<<< HEAD
    s_in(i).x = s_in(i).x(:);
    [s_in(i).x, ind] = sort(s_in(i).x);
    s_in(i).y = s_in(i).y(:);
    s_in(i).y = s_in(i).y(ind);
    s_in(i).e = s_in(i).e(:);
    s_in(i).e = s_in(i).e(ind);
    if ~isempty(s_in(i).yfit)
        s_in(i).yfit = s_in(i).yfit(:);
        s_in(i).yfit = s_in(i).yfit(ind);
=======
    for j = 1:length(f)
        if isfield(s_in(i),f{j})
            s(i).(f{j}) = s_in(i).(f{j});
            if any(strcmpi(f{j},{'x_label','y_label','datafile'}))
                if isnumeric(s(i).(f{j}))
                    s(i).(f{j}) = num2str(s(i).(f{j}));
                end
            end
        end
>>>>>>> 636511af990e58b16bd962036363f5ae877ec4b8
    end
    validate(s(i));
end

s_out = reshape(s,sz);
