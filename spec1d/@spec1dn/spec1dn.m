classdef spec1dn < spec1d
    %SPECRES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=protected)
        resolution
        polarisation
    end
    
    properties (Hidden=true,Dependent=true,SetAccess=protected)
        crystal
    end
    
    methods
        function obj = spec1dn(varargin)

            if nargin == 0
                return;
            end
            
            if nargin ~= 0
                if isa(varargin{1},'spec1d')
                    s = varargin{1};
                    obj = obj.copy(s);
                    if length(varargin) == 2
                        switch class(varargin{2})
                            case 'char'
                                if exist(varargin{2},'file')
                                    obj.datafile = varargin{2};
                                end
                            case 'struct'
                                if isfield(varargin{2},'EXP')
                                    obj.resolution = varargin{2};
                                end
                            otherwise
                                warning('Could not assign datafile!')
                        end
                    end
                end
                if any(nargin == [5 6 7])
                    loader = varargin{1};
                    filename = varargin{2};
                    X = varargin{3};
                    Y = varargin{4};
                    M = varargin{5};
                    if nargin >= 6
                        if ischar(varargin{6})
                            if ~isempty(strfind(varargin{6},'PAL'))
                                opt = varargin{6};
                                j = 1;
                                cont = true;
                                while cont;
                                    try
                                        opt_F = sprintf('%s,F=%i',opt,j);
                                        temp(j) = loads(loader,sprintf('%s,X=%s,Y=%s,M=%s,%s',filename,X,Y,M,opt_F));
                                        j = j +1;
                                    catch
                                        cont = false;
                                    end
                                end
                            else
                                temp = loads(loader,sprintf('%s,X=%s,Y=%s,M=%s,%s',filename,X,Y,M,varargin{6}));
                            end
                            if nargin == 7
                                if ~exist('temp','var')
                                    error('Could not load file')
                                end
                                temp = temp*varargin{7};
                            end
                        else
                            temp = loads(loader,sprintf('%s,X=%s,Y=%s,M=%s',filename,X,Y,M))*varargin{6};
                        end
                    else
                        temp = loads(loader,sprintf('%s,X=%s,Y=%s,M=%s',filename,X,Y,M));
                    end
                    for i = 1:length(temp);
                        obj(i) = spec1dn();
                        obj(i) = obj(i).copy(temp);
                        if isempty(obj(i).datafile)
                            [~, temp2] = fileparts(filename);
                        else
                            [~, temp2] = fileparts(obj(i).datafile);
                        end
                        obj(i).datafile = fullfile(fileparts(filename),temp2);
                    end
                end
            end
            for i = 1:length(obj)
                if ~isempty(obj(i).datafile) && isempty(obj(i).resolution)
                    obj(i).resolution = ResLibCal('silent','open',obj(i).datafile);
                    obj(i).resolution.dqh = [obj(i).resolution.ResCal.DQH,...
                        obj(i).resolution.ResCal.DQK,...
                        obj(i).resolution.ResCal.DQL,...
                        obj(i).resolution.ResCal.DEN];
                end
                if exist('cont','var')
                    obj(i) = obj(i).addPolarisation();
                end
            end
        end
        
        function value = get.crystal(obj)
            if ~isempty(obj.resolution)
                obj = madrlp(obj);
                value = obj.crystal;
            else
                error('Resoltion has not been set. Use obj.setResolution(cfg)')
            end
        end
        
        function disp(obj)
            disp@spec1d(obj)
            if length(obj) == 1
                if isempty(obj.resolution)
                    fprintf('Resoltion has not been set. Use obj.setResolution(cfg).\n')
                else
                    fprintf('Resoltion computed for Qh=%0.03f, Qk=%0.03f, Ql=%0.03f, En=%0.03f\n',...
                        obj.resolution.EXP.QH,obj.resolution.EXP.QK,obj.resolution.EXP.QL,obj.resolution.EXP.W)
                end
                if isempty(obj.polarisation)
                    fprintf('This is unpolarised data.\n')
                else
                    fprintf('The data is polarised, but consisting of %i channel.\n',length(obj.polarisation))
                end
            else
                if isempty(obj(1).resolution)
                    fprintf('Resoltion has not been set. Use obj.setResolution(cfg).\n')
                else
                    fprintf('Resoltion computed for Qh=%0.03f-%0.03f, Qk=%0.03f-%0.03f, Ql=%0.03f-%0.03f, En=%0.03f-%0.03f\n',...
                        min(arrayfun(@(x)x.EXP.QH,[obj.resolution])),max(arrayfun(@(x)x.EXP.QH,[obj.resolution])),...
                        min(arrayfun(@(x)x.EXP.QK,[obj.resolution])),max(arrayfun(@(x)x.EXP.QK,[obj.resolution])),...
                        min(arrayfun(@(x)x.EXP.QL,[obj.resolution])),max(arrayfun(@(x)x.EXP.QL,[obj.resolution])),...
                        min(arrayfun(@(x)x.EXP.W ,[obj.resolution])),max(arrayfun(@(x)x.EXP.W ,[obj.resolution])))
                end
                if isempty(obj(1).polarisation)
                    fprintf('This is unpolarised data.\n')
                else
                    fprintf('The data is polarised with %i channels.\n',sum(arrayfun(@(x)~isempty(x.polarisation),obj)))
                end
            end
        end
        
        function obj = addResolution(obj,cfg)
            if isstruct(cfg)
                if isfield(cfg,'EXP')
                    obj.resolution = cfg;
                end
                return
            end
            if ischar(cfg)
                if exist(cfg,'file')
                    obj.resolution = ResLibCal('silent','open',cfg);
                    obj = madrlp(obj);
                    obj.resolution.dqh = [obj.resolution.ResCal.DQH,...
                        obj.resolution.ResCal.DQK,...
                        obj.resolution.ResCal.DQL,...
                        obj.resolution.ResCal.DEN];
                else
                    error('File does not exist')
                end
            else
                error('Unknown option')
            end
        end
        
        function obj = addPolarisation(obj,varargin)
            obj.polarisation = obj.createPolarisation(varargin{:});
            obj.userdata.addon = struct('polarised',true,'channel',obj.polarisation.channel);
        end
        
    end
    
    methods (Access = protected)
        function resolutionMouseButtonPress(obj,obj_a,event)
            
            obj = obj(obj.findID(obj_a.UserData.ID));
                        
            if isempty(obj.resolution)
                return
            else
                [~, ind] = min(abs(obj.x - event.IntersectionPoint(1)));
                if abs(obj.y(ind)-event.IntersectionPoint(2)) > 0.05*diff(get(get(obj_a,'Parent'),'YLim'))
                    return
                end
                
                if obj.resolution.dqh(1)
                    h = obj.x(ind);
                else
                    h = obj.resolution.EXP.QH;
                end
                if obj.resolution.dqh(2)
                    k = obj.x(ind);
                else
                    k = obj.resolution.EXP.QK;
                end
                if obj.resolution.dqh(3)
                    l = obj.x(ind);
                else
                    l = obj.resolution.EXP.QL;
                end
                if obj.resolution.dqh(4)
                    e = obj.x(ind);
                else
                    e = obj.resolution.EXP.W;
                end
                obj = obj.setQ(h,k,l,e);
                temp = ResLibCal('silent','view2',obj.resolution);
            end
        end
        
        function f = createPolarisation(obj,varargin)
            
            p = inputParser;
            p.addOptional('channel',[])
            p.addOptional('F1',0)
            p.addOptional('F2',0)
            
            p.parse(varargin{:})
            
            f = p.Results;
        end
        
    end
end
