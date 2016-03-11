function out = fieldnames_s(obj)

out = arrayfun(@(x)arrayfun(@(y)y.Name,x.PropertyList,'UniformOutput',0),eval(sprintf('?%s',class(obj))),'UniformOutput',0);
out = out{:};
end