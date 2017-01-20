function varargout = fits(obj,varargin)
%FITS Summary of this function goes here
%   Detailed explanation goes here

test = strcmpi('noRes',varargin);
if any(test)
    varargin(test) = [];
else
    hkl = cell(4,1);
    for i = 1:length(obj)
        cfg = ResCalFit_ini(obj(i));
        tmp = cell(4,1);
        if isfield(cfg,'SQW') || isfield(cfg,'BG')
            if isempty(cfg.SQW) || isempty(cfg.BG)
                error('SQW and Bg can not be empty')
            else
                if cfg.dqh(1)
                    tmp{1} = obj(i).x;
                else
                    tmp{1} = cfg.EXP.QH;
                end
                if cfg.dqh(2)
                    tmp{2} = obj(i).x;
                else
                    tmp{2} = cfg.EXP.QK;
                end
                if cfg.dqh(3)
                    tmp{3} = obj(i).x;
                else
                    tmp{3} = cfg.EXP.QL;
                end
                if cfg.dqh(4)
                    tmp{4} = obj(i).x;
                else
                    tmp{4} = cfg.EXP.W;
                end
                m = max(cellfun(@length,tmp));
                if length(tmp{1}) == 1
                    tmp{1} = repmat(tmp{1},1,m);
                end
                if length(tmp{2}) == 1
                    tmp{2} = repmat(tmp{2},1,m);
                end
                if length(tmp{3}) == 1
                    tmp{3} = repmat(tmp{3},1,m);
                end
                if length(tmp{4}) == 1
                    tmp{4} = repmat(tmp{4},1,m);
                end
            end
        else
            error('Use addSQW')
        end
        hkl{1} = [hkl{1}(:); tmp{1}(:)];
        hkl{2} = [hkl{2}(:); tmp{2}(:)];
        hkl{3} = [hkl{3}(:); tmp{3}(:)];
        hkl{4} = [hkl{4}(:); tmp{4}(:)];
    end
    ResCalFit_ini(cfg,hkl{1},hkl{2},hkl{3},hkl{4})
end

varargout = cell(nargout,1);
[varargout{:}] = fits@spec1d(obj,varargin{:});

end

