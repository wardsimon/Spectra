function obj = matparser(obj, varargin)

vals = obj.p_vals;

inpForm.fname  = {'param' 'mat' 'vals'};
inpForm.defval = {[]      []     []     };
inpForm.size   = {[1 -1]  [1 -1] []     };

if nargin >= 2 || isempty(obj.parser)
    out = parse_it(inpForm,varargin{:});
else
    for i = 1:length(inpForm.fname)
        out.(inpForm.fname{i}) = inpForm.defval{i};
    end
    out.mat = obj.parser;
    out.vals = varargin{:};
end


if ~isempty(out.param)
    vals = out.param;
elseif ~ isempty(out.mat) && isstruct(vals)
    if length(out.mat) == length(out.vals)
        for i = 1:length(out.mat)
            if isfield(vals,out.mat{i})
                vals.(out.mat{i}) = out.vals{i};
            else
                try
                    eval(sprintf('vals.%s = out.vals{%i}',out.mat{i},i))
                catch
                    warning('Can not set it')
                end
            end
        end
    end
end

obj.p_vals = vals;

end

function out = parse_it(format,varargin)

for i = 1:length(format.fname)
    raw.(format.fname{i}) = format.defval{i};
end

if (nargin>2) && (mod(nargin,2) == 1)
    nPar = nargin-1;
    for ii = 1:2:nPar
        raw.(varargin{ii}) = varargin{ii+1};
    end
elseif nargin == 2
    raw = varargin{1};
else
    error('sw:sw_readparam:WrongNumberOfInput','Wrong number of input options!');
end

out = mergestruct(format,raw);

end