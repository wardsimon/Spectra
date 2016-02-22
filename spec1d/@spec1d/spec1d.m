function s = spec1d(a,varargin)
%
% function spec1d(varargin)
%
% SPEC1D/spec1d Create a spec1d spectra object
%
% Usage:    1. s = spec1d(x,y,e);
%           3. s = spec1d(x,y,,e,'datafile',file,'x_label',xlabel,'y_label',ylabel);
%
% Update : Start GPU work.
% Simon Ward 25/01/2016

    p = inputParser;
    
    if nargin == 0
        p.addParamValue('x',[],@(x) isnumeric(x) & isreal(x))
        p.addParamValue('y',[],@(x) isnumeric(x) & isreal(x))
        p.addParamValue('e',[],@(x) isnumeric(x) & isreal(x) & all(x >= 0))
        p.addParamValue('yfit',[],@(x) isnumeric(x) & isreal(x))
        p.addParamValue('x_label',[],@(x) ischar(x) | isempty(x))
        p.addParamValue('y_label',[],@(x) ischar(x) | isempty(x))
        p.addParamValue('datafile',[],@(x) ischar(x) | isempty(x))
        p.parse();
        s = class(p.Results,'spec1d');
        return
    end
    
    if isstruct(a) || isa(a,'spec1d')
        if isa(a,'spec1d')
            a = struct(a);
        end
        p.addParamValue('x',[],@(x) isnumeric(x) & isreal(x))
        p.addParamValue('y',[],@(x) isnumeric(x) & isreal(x))
        p.addParamValue('e',[],@(x) isnumeric(x) & isreal(x) & all(x >= 0))
        p.addParamValue('yfit',[],@(x) isnumeric(x) & isreal(x))
        p.addParamValue('x_label',[],@(x) ischar(x) | isempty(x))
        p.addParamValue('y_label',[],@(x) ischar(x) | isempty(x))
        p.addParamValue('datafile',[],@(x) ischar(x) | isempty(x))
    else
        if isempty(varargin)
            varargin{1} = zeros(size(a));
        end
        if length(varargin) == 1
            varargin{2} = zeros(size(a));
        end
        p.addRequired('x',@(x) isnumeric(x) & isreal(x))
        p.addRequired('y',@(x) isnumeric(x) & isreal(x))
        p.addRequired('e',@(x) isnumeric(x) & isreal(x) & all(x >= 0))
        p.addOptional('yfit',[],@(x) isnumeric(x) & isreal(x))
        p.addOptional('x_label',[],@ischar)
        p.addOptional('y_label',[],@ischar)
        p.addOptional('datafile',[],@ischar)
    end
    
    if length(a) > 1 && isstruct(a)
        si = size(a);
        a = a(:);
        for i = 1:length(a)
            s_in(i) = spec1d(a(i));
        end
        s = reshape(s_in,si);
        return
    else
        p.parse(a,varargin{:});
        s_in = p.Results;
    end
    
    for i = 1:length(s_in)
        s_in(i).x = s_in(i).x(:);
        s_in(i).y = s_in(i).y(:);
        if length(s_in(i).e(:))==1
            s_in(i).e = s_in(i).e(:)*ones(size(s_in(i).x(:)));
        end
        s_in(i).e = s_in(i).e(:);
        s_in(i).yfit = s_in(i).yfit(:);
        
        if any(diff([length(s_in(i).x) length(s_in(i).y) length(s_in(i).e)]))
            warning('x, y and e are not of the same length!\nLength x:\t%i\nLength y:\t%i\nLength e:\t%i',length(s_in(i).x),length(s_in(i).y),length(s_in(i).e))
        end
        
        if (length(s_in(i).x) > 1E5) && sdext.getpref('gpuArray').val && checkGpuMemory(s_in(i).x)
           s_in(i).x = gpuArray(s_in(i).x); 
           s_in(i).y = gpuArray(s_in(i).y);
           s_in(i).e = gpuArray(s_in(i).e);
        else
           if isa(s_in(i).x,'gpuArray') % One of the conditions is not met, pull back
               s_in(i).x = gather(s_in(i).x);
               s_in(i).y = gather(s_in(i).y);
               s_in(i).e = gather(s_in(i).e);
           end
        end
        s(i) = class(s_in(i),'spec1d');
    end 
end

function OK = checkGpuMemory(a)
    d = gpuDevice;
    memA = d.AvailableMemory;  % Check available memory
    memD = 3*length(a)*8; % We will use this much memory
    if memD < memA
        OK = 1;
    else
        OK = 0;
    end
end