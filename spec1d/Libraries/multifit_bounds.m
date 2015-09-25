function bounds_out = multifit_bounds(flag,bounds)
%% Multifit bounds procedure
% flag = 0, 1, 2, 3
% bounds = cell of 

%----- dpin now has 3 possible states: 0=fixed, 1=the same for all spec1d
%      objects, 2=fitted independently in all spec1d objects. I have now to
%      create the actual pin and dp vectors

%----- Also, param_keep will store the order of the parameters
global  x_per_spec param_keep

bounds_out=zeros((length(flag)+length(find(flag>1))*(length(x_per_spec)-1)),2);

ll = 1;
for j = 1:length(flag)
    if flag(j) == 2
        if size(bounds{j},2) == 2
            bounds{j} = repmat(bounds{j},length(ll:(ll+length(x_per_spec)-1)),1);
        end
        bounds_out(ll:(ll+length(x_per_spec)-1),:) = bounds{j};
        ll = ll+length(x_per_spec);
    elseif flag(j) == 3
        if size(bounds{j},2) == 2
            bounds{j} = repmat(bounds{j},length(ll:(ll+length(x_per_spec)-1)),1);
        end
        bounds_out(ll:(ll+length(x_per_spec)-1),:) = bounds{j};
        ll = ll+length(x_per_spec);
    elseif (flag(j) == 1) || (flag(j) == 0)
        bounds_out(ll,:) = bounds{j};
        ll = ll+1;
    else
        disp('Please use following notation: 0 = fixed, 1 = same for all objects, 2 = independently fitted in all objects, 3 = fixed and different for all objects')
    end
end
