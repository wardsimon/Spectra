function s_out = loadobj(s_in)
%% This function fixes the conversion between the original spec1d and the post 2013 version.
sz = size(s_in);

s_in = s_in(:);

s = spec1d();
f = fieldnames(s);
s = repmat(s,size(s_in));

for i = 1:length(s_in)
    for j = 1:length(f)
        if isfield(s_in(i),f{j})
            s(i).(f{j}) = s_in(i).(f{j});
            if any(strcmpi(f{j},{'x_label','y_label','datafile'}))
                if isnumeric(s(i).(f{j}))
                    s(i).(f{j}) = num2str(s(i).(f{j}));
                end
            end
        end
    end
    validate(s(i));
end

s_out = reshape(s,sz);
