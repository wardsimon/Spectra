function s_out = loadobj(s_in)
%% This function fixes the conversion between the original spec1d and the simon version.


s = size(s_in);

s_in = s_in(:);

for i = 1:length(s_in)
    s_out(i) = spec1d(s_in(i));
end

s_out = reshape(s_out,s);
