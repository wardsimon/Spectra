function rPref = getpref(prefName)
% returns specnd global preferences
%
% rPref = ndext.getpref
%
% Returns the names, values and labels of each preferences. Default values
% are returned, where no value is saved. rPref is a struct with field names
% 'name', 'label' and 'val'. Each field is a cell.
%
% rPref = ndext.getpref(pName)
%
% Returns only the requested specnd preference name, value and label. Each
% field contains the requested value.
%
% rPref = ndext.getpref('default')
%
% Returns the default names, values and labels of each preferences.
%
% See also ndext.setpref.
%

% binmethod	center/weight/sliding	Binning method
% autosort	true/false	Automatically sort data and coordinates after certain operations
% ach	[], 1, 2, 3, ...	Index of the active channel
% importfun	function_handle	Default import function
% datapath	string	Path to the data folder for import
% filename	string	File name format to fast file import
% emptyval	NaN,-inf,inf,0,...	Return value for empty data point
% mon	double	Global monitor

% default values
dn = {      'fid'       'binmethod' 'bintype'   'autosort'  'ach'   };
dv = {      1           'center'    'center'    'true'      []      };
dn = [dn {  'importfun' 'datapath'  'filename'  'emptyval'  'mon'   }];
dv = [dv {  []          ''          ''          nan         1       }];
dn = [dn {     'errfun'    }];
dv = [dv {     @sqrt       }];

dl = {...
    'file identifier for text output'...
    'signal binning method: center/weight/sliding'...
    'type of bin in histogram mode: center/boundary'...
    'automatically sort coordinates when necessary'...
    'index of active channel'...
    'handle of the data import function'...
    'path to the data folder'...
    'file name format to fast file import'...
    'signal value for empty pixels'...
    'global monitor value'...
    'function that calculates the error bar from a given signal value'...
    };

dPref = struct('val',{},'name',{},'label',{});

[dPref(1:numel(dv),1).name] = dn{:};
[dPref(:).label]          = dl{:};
[dPref(:).val]            = dv{:};

% get stored preferences
sPref = getpref('specnd_global');

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
        error('specnd:getndpref:WrongName','The requested specnd preference does not exists!');
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