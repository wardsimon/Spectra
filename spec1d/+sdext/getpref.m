function rPref = getpref(prefName)
% returns spec1d global preferences
%
% rPref = sdext.getpref
%
% Returns the names, values and labels of each preferences. Default values
% are returned, where no value is saved. rPref is a struct with field names
% 'name', 'label' and 'val'. Each field is a cell.
%
% rPref = sdext.getpref(pName)
%
% Returns only the requested specnd preference name, value and label. Each
% field contains the requested value.
%
% rPref = sdext.getpref('default')
%
% Returns the default names, values and labels of each preferences.
%
% See also sdext.setpref.
%
% Branched from specnd 
%
% default values

path_ind = strfind(mfilename('fullpath'),'/');
path_d =  mfilename('fullpath');
path_d = path_d(1:path_ind(end-3));

dn = { 'libroot' 'experimental' 'doLog','gpuArray'};
dv = { path_d    0              2        0};
dl = {...
    'Where Spectra files are stored',...
    'Enable experiemntal features',...
    'Log all interactions',...
    'Should we run on the GPU'};

dPref = struct('val',{},'name',{},'label',{});

[dPref(1:numel(dv),1).name] = dn{:};
[dPref(:).label]          = dl{:};
[dPref(:).val]            = dv{:};

% get stored preferences
sPref = getpref('mtools');

if (nargin>0)
    if strcmp(prefName,'default')
        % return default preference values
        rPref = dPref;
        return
    end
    
    % if a specific value is requested, check if it exist in the default value
    % list
    
    iPref = find(strcmp(prefName,{dPref(:).name}),1);
    if isempty(iPref)
        error('spectra:getsdpref:WrongName','The requested spectra preference does not exist!');
    end
    
    % if a specific value is requested and it exists, return it
    rPref = dPref(iPref);
    
    if isfield(sPref,prefName)
        rPref.val = sPref.(prefName);
    end
    
    return
else
    % return all stored values
    rPref = dPref;
    % overwrite default values for existing preferences
    if ~isempty(sPref)
        fPref = fieldnames(sPref);
        for ii = 1:numel(fPref)
            rPref(strcmp(fPref{ii},{dPref(:).name})).val = sPref.(fPref{ii});
        end
    end
    
end

end