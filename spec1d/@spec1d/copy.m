function out = copy(obj,objnew)

if nargin == 1
    objnew = obj;
end

fn1 = fieldnames_s(obj);
if isstruct(objnew)
    fn2 = fieldnames(objnew);
    super = 0;
else
    fn2 = fieldnames_s(objnew);
    out = arrayfun(@(y) y.SuperclassList,eval(sprintf('?%s',class(objnew))),'UniformOutput',0);
    super = ~isempty(out{:}); % Determine if this is a spec1d or inherited class.
end

if super
    out = feval(class(objnew));
else
    out = feval(class(obj));
end

for i = 1:length(fn2)
    if any(strcmp(fn2{i},fn1))
        try
            out.(fn2{i}) = objnew.(fn2{i});
        catch
            % We can not set the unique identifier field.
        end
    end
end

out.fitdata.specID = out.ident;

