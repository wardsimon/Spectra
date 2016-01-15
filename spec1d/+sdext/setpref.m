function setpref(prefName, value)
% sets spec1d global preferences
%
% sdext.setpref(prefName, value)
%
% Sets the value of the prefName specnd global preferences.
%
% sdext.setpref('default')
%
% Resets all preferences values to the default one.
%
% See also sdext.getpref.
%
% Branched from specnd
%

if strcmp(prefName,'default')
    if ispref('mtools')
        rmpref('mtools');
    end
    return
end

% check if the preference name exists
dPref = sdext.getpref('default');

iPref = find(strcmp(prefName,{dPref(:).name}),1,'first');
if ~isempty(iPref)
    % check if the preferences label contains a choice string
    str0 = strsplit(dPref(iPref).label,' ');
    opt0 = strsplit(str0{end},'/');
    
    if numel(opt0) > 1
        % there is a choice of different string options
        if ~ischar(value) || ~any(strcmp(value,opt0))
            error('spectra:setpref:WrongInput',['The selected preference has a restricted choice: ' str0{end} '!'])
        end
        setpref('mtools',prefName,value);
    else
        % the value has to be a scalar
        % TODO check for other type of values
        setpref('mtools',prefName,value);
    end
    
else
    error('spectra:setpref:WrongName','The given name is not a valid spectra global preferences!');
end

end