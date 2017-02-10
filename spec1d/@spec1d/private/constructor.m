function s = constructor(a,varargin)
% Parses all the arguments for a spec1d object
%
% s_out = SPEC1D(varargin)
%
% Input:
% 
% x,y,e     Vectors of real numbers where e is positive and non-zero. e can
%           be a single value, replicates for all datapoints.
% s_in      Structure with fields x, y, e or an existing spec1d object.
%
% Optional:
% 
% fitdata   Fit data stored in the @specfit class.
% yfit      y-points of a fit which have been evaluated at the x-points of
%           s_in.
% x_label   Text field corresponding to a label which goes on the x-axis of
%           a plot.
% y_label   Text field corresponding to a label which goes on the y-axis of
%           a plot.
% datafile  Text field corrsponding to where the data was obtained. Will
%           also be used as a title for plotting.
% userdata  A structure which stores additional data. This can be set as
%           required and has no set fields.
%
% This function creates a spec1d object from inputs given. 
%
% Example:    
%   s_out = spec1d(x,y,e);
%   s_out = spec1d(x,y,e,'datafile',file,'x_label',xlabel,'y_label',ylabel);
%

p = inputParser;

if nargin == 0
    p.addParamValue('x',[],@(x) validateattributes(x,{'numeric','gpuArray'},{'vector','real','finite','nonnan'},mfilename,'x'))
    p.addParamValue('y',[],@(x) validateattributes(x,{'numeric','gpuArray'},{'vector','real','finite','nonnan'},mfilename,'y'))
    p.addParamValue('e',[],@(x) validateattributes(x,{'numeric','gpuArray'},{'nonnegative','vector','real','finite','nonnan'},mfilename,'e'))
    p.addParamValue('fitdata',specfit(),@(x)isa(x,'specfit'))
    p.addParamValue('yfit',[],@(x) validateattributes(x,{'numeric','gpuArray'},{'real','finite','nonnan','2d'},mfilename,'yfit'))
    p.addParamValue('x_label','',@(x) isempty(x) || validateattributes(x,{'char'},{},mfilename,'x_label'))
    p.addParamValue('y_label','',@(x) isempty(x) || validateattributes(x,{'char'},{},mfilename,'y_label'))
    p.addParamValue('datafile','',@(x) isempty(x) || validateattributes(x,{'char'},{},mfilename,'datafile'))
    p.addParamValue('userdata',[])
    p.parse();
    s = p.Results;
    return
end

if isstruct(a) || isa(a,'spec1d')
    if isa(a,'spec1d')
        temp = a;
        a = struct;
        pr = properties(temp);
        for i = 1:length(pr)
            a.(pr{i}) = temp.(pr{i});
        end
    end
    p.addParamValue('x',[],@(x) validateattributes(x,{'numeric','gpuArray'},{'vector','real','finite','nonnan'},mfilename,'x'))
    p.addParamValue('y',[],@(x) validateattributes(x,{'numeric','gpuArray'},{'vector','real','finite','nonnan'},mfilename,'y'))
    p.addParamValue('e',[],@(x) validateattributes(x,{'numeric','gpuArray'},{'nonnegative','vector','real','finite','nonnan'},mfilename,'e'))
    p.addParamValue('yfit',[],@(x) validateattributes(x,{'numeric','gpuArray'},{'real','finite','nonnan','2d'},mfilename,'yfit'))
    p.addParamValue('fitdata',specfit(),@(x)isa(x,'specfit'))
    p.addParamValue('x_label','',@(x) validateattributes(x,{'char'},{'2d'},mfilename,'x_label'))
    p.addParamValue('y_label','',@(x) validateattributes(x,{'char'},{'2d'},mfilename,'y_label'))
    p.addParamValue('datafile','',@(x) validateattributes(x,{'char'},{'2d'},mfilename,'datafile'))
    p.addParamValue('userdata',[])
else
    if isempty(varargin)
        varargin{1} = zeros(size(a));
    end
    if length(varargin) == 1
        varargin{2} = zeros(size(a));
    end
    p.addRequired('x',@(x) validateattributes(x,{'numeric','gpuArray'},{'vector','real','finite','nonnan'},mfilename,'x'))
    p.addRequired('y',@(x) validateattributes(x,{'numeric','gpuArray'},{'vector','real','finite','nonnan'},mfilename,'y'))
    p.addRequired('e',@(x) validateattributes(x,{'numeric','gpuArray'},{'nonnegative','vector','real','finite','nonnan'},mfilename,'e'))
    p.addParamValue('yfit',[],@(x) validateattributes(x,{'numeric','gpuArray'},{'real','finite','nonnan','2d'},mfilename,'yfit'))
    p.addParamValue('fitdata',specfit(),@(x)isa(x,'specfit'))
    p.addParamValue('x_label','',@(x)validateattributes(x,{'char'},{'2d'},mfilename,'x_label'))
    p.addParamValue('y_label','',@(x)validateattributes(x,{'char'},{'2d'},mfilename,'y_label'))
    p.addParamValue('datafile','',@(x)validateattributes(x,{'char'},{'2d'},mfilename,'datafile'))
    p.addParamValue('userdata',[])
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
    [s_in(i).x, ind]= sort(s_in(i).x(:));
%     s_in(i).x = s_in(i).x(:);
    s_in(i).y = reshape(s_in(i).y(ind),[],1);
%     s_in(i).y = s_in(i).y(:);
    if length(s_in(i).e(:)) == 1
        s_in(i).e = s_in(i).e(:)*ones(size(s_in(i).x(:)));
    end
    s_in(i).e = reshape(s_in(i).e(ind),[],1);
%     s_in(i).e = s_in(i).e(:);
    if ~isempty(s_in(i).yfit)
        s_in(i).yfit = reshape(s_in(i).yfit(ind),[],1);
%         s_in(i).yfit = s_in(i).yfit(:);
    else
        s_in(i).yfit = s_in(i).yfit(:);
    end
    if isfield(s_in(i).userdata,'rind')
        % Dont f*ck around with indexing!
        if isempty(s_in(i).userdata.rind)
            [~, s_in(i).userdata.rind] = sort(ind);
        % Check if this is multifit call or something.    
        elseif length(s_in(i).userdata.rind) ~= length(s_in(i).x)
            [~, s_in(i).userdata.rind] = sort(ind);
        end
    else
        [~, s_in(i).userdata.rind] = sort(ind);
    end
    
    if any(diff([length(s_in(i).x) length(s_in(i).y) length(s_in(i).e)]))
        warning('x, y and e are not of the same length!\nLength x:\t%i\nLength y:\t%i\nLength e:\t%i',length(s_in(i).x),length(s_in(i).y),length(s_in(i).e))
    end
    
    if (length(s_in(i).x) > sdext.getpref('minGPULength').val) && sdext.getpref('gpuArray').val && checkGpuMemory(s_in(i))
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
end
s = s_in;
end