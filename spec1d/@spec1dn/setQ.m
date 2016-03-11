function obj = setQ(obj,varargin)

p = inputParser;
p.addOptional('H',obj.resolution.EXP.QH,@(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','2d'},mfilename,'H'))
p.addOptional('K',obj.resolution.EXP.QK,@(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','2d'},mfilename,'K'))
p.addOptional('L',obj.resolution.EXP.QL,@(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','2d'},mfilename,'L'))
p.addOptional('E',obj.resolution.EXP.W,@(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','2d'},mfilename,'E'))
if isfield(obj.resolution,'dqh')
    p.addOptional('scan',obj.resolution.dqh,@(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','2d'},mfilename,'scan'))
else
    dqh = [obj.resolution.ResCal.DQH obj.resolution.ResCal.DQK obj.resolution.ResCal.DQL obj.resolution.ResCal.DEN];
    p.addOptional('scan',dqh,@(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','2d'},mfilename,'scan'))
end

p.parse(varargin{:});
results = p.Results;

obj.resolution.EXP.QH = results.H; obj.resolution.QH = results.H; 
obj.resolution.EXP.QK = results.K; obj.resolution.QK = results.K; 
obj.resolution.EXP.QL = results.L; obj.resolution.QL = results.L; 
obj.resolution.EXP.W = results.E;  obj.resolution.W = results.E; 

if results.scan(1)
    obj.resolution.ResCal.DQH = results.scan(1);
end
if results.scan(2)
    obj.resolution.ResCal.DQK = results.scan(2);
end
if results.scan(3)
    obj.resolution.ResCal.DQL = results.scan(3);
end
if results.scan(4)
    obj.resolution.ResCal.DEN = results.scan(4);
end
obj.resolution.dqh = results.scan;


temp = ResLibCal(obj.resolution.EXP,'silent','compute');
obj.resolution = mergestruct(obj.resolution,temp);
end