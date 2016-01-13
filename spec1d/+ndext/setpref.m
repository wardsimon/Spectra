function setpref(prefName, value)
% sets specnd global preferences
%
% ndext.setpref(prefName, value)
%
% Sets the value of the prefName specnd global preferences.
%
% ndext.setpref('default')
%
% Resets all preferences values to the default one.
%
% See also ndext.getpref.
%

if strcmp(prefName,'default')
    if ispref('specnd_global')
        rmpref('specnd_global');
    end
    return
end

% check if the preference name exists
dPref = ndext.getpref('default');

iPref = find(strcmp(prefName,{dPref(:).name}),1,'first');
if ~isempty(iPref)
    % check if the preferences label contains a choice string
    str0 = strsplit(dPref(iPref).label,' ');
    opt0 = strsplit(str0{end},'/');
    
    if numel(opt0) > 1
        % there is a choice of different string options
        if ~ischar(value) || ~any(strcmp(value,opt0))
            error('ndext:setpref:WrongInput',['The selected preference has a restricted choice: ' str0{end} '!'])
        end
        setpref('specnd_global',prefName,value);
    else
        % the value has to be a scalar
        % TODO check for other type of values
        setpref('specnd_global',prefName,value);
    end
    
else
    error('ndext:setpref:WrongName','The given name is not a valid specnd global preferences!');
end

end