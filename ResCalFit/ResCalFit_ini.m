function varargout = ResCalFit_ini(file,varargin)
%% ResCalFit_ini function scan = ResCalFit_ini(file,varargin)
% Initialisation script for fitting triple-axis data using ResLibCal as a back-end.
% !!! This function typically needs 2 calls !!!
%
% First call
% file can be a ResLibCal file, scan file, trixfit file, ResLib file etc...
% varargin is empty.
%
% Get the initial configuration by just supplying the file.
%
% Second Call
% Supply the configuration obtained from the first call as the first argument.
%
% Arguments (Required) after file call....
% file = configuration from previous call.
% varargin{1} = Qh;  Single value or vector
% varargin{2} = Qk;  Single value or vector
% varargin{3} = Ql;  Single value or vector
% varargin{4} = E;   Single value or vector
% varargin{5} = Cross-Section file;  'Text'
% varargin{6} = Background file;     'Text'
%
% S. Ward (simon.ward@psi.ch), September 2015
% Version 1.3 (18.01.2016) - Support for GPU

global ResFitScn

if isempty(file)
    cfg = ResLibCal('compute');
else
    if ischar(file)
        cfg = ResLibCal('open',file);
    elseif isstruct(file)
        if isfield(file,'EXP')
            cfg = ResLibCal(file,'compute');
        end
    elseif isa(file,'specres')
        cfg = file.resolution;
    end
end

% Process arguments
if isempty(varargin)
    varargout{1} = cfg;
    return
else
    if length(varargin) < 4
        error('You must supply arguments H,K,L and E');
    end
    % Sort out H,K,L,E of various lengths
    larg = cellfun(@length,varargin(1:4));
    Q = zeros(max(larg),4);
    for i = 1:4
        if length(varargin{i})==1
            Q(:,i) = varargin{i};
        else
            if length(varargin{i}) ~= max(larg)
                error('Wrong length');
            else
                Q(:,i) = varargin{i};
            end
        end
    end
    % Add the cross-section and background file.
    if length(varargin) == 4
        if isfield(cfg,'SQW') && isfield(cfg,'BG')
            if isempty(cfg.SQW) && isempty(cfg.BG)
                error()
            end
        else
            error()
        end
    else
        if length(varargin) <= 7
            cfg.SQW = varargin{5};
            cfg.BG = varargin{6};
            if length(varargin) ==7
                cfg.EXP.NMC = varargin{7};
            end
        end
    end
end

% Create the scan
if length(cfg) == 1
    ResFitScn = repmat(cfg,1,max(larg));
else
    ResFitScn = cfg;
    if nargout == 1
        varargout{1} = ResFitScn;%arrayfun(@(x) ResLibCal(x,'compute'),scan);
    end
    return
end

% Fill in the scan with the correct Q points and re-compute
for i = 1:length(ResFitScn)
    ResFitScn(i).EXP.QH = Q(i,1); ResFitScn(i).ResCal.QH = Q(i,1);
    ResFitScn(i).EXP.QK = Q(i,2); ResFitScn(i).ResCal.QK = Q(i,2);
    ResFitScn(i).EXP.QL = Q(i,3); ResFitScn(i).ResCal.QL = Q(i,3);
    ResFitScn(i).EXP.W  = Q(i,4); ResFitScn(i).ResCal.EN = Q(i,4);
    ResFitScn(i).resolution.HKLE = Q(i,:);
    ResFitScn(i) = ResLibCal(ResFitScn(i),'silent','compute');
end

% ResFitScn = arrayfun(@updateMC,scan);

if nargout == 1
    varargout{1} = ResFitScn;%arrayfun(@(x) ResLibCal(x,'compute'),scan);
end
end