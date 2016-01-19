function s_in = updateMC(s_in)
% ResLibCal_RM2clouds: creates an axis system of Monte-Carlo points
%   which represent the resolution function in [abc] and [xyz] frames.
%
%   EXP: structure containing application configuration. When not given, extracts
%        it from the main window.
%   resolution: the resolution computation structure, for both [abc] and [xyz] frames
%
% Returns:
%   resolution.rlu.cloud: cloud of points 4D axes as { H,K,L,W } in [ABC] frame
%   resolution.spec.cloud: cloud of points 4D axes as { H,K,L,W } in [xyz] frame

% insprired from ResCal5/mc_conv
% this routine uses: in [abc|xyz] frame: RM (that is along [ABC|xyz] axes)
%                    the HKLE location, converted into [ABC|xyz] frame with hkl2Frame)

% accuracy obtained from McStas estimates:
% NMC=1000  -> 10%
% NMC=10000 -> 2.5%

try
    % This assumes the function was called from ResFitCal
    runGPU = evalin('caller','runGPU');
catch
    runGPU = sdext.getpref('gpuArray').val;
end
    
if isfield(s_in,'EXP') && isfield(s_in,'resolution')
    resolution = s_in.resolution;
    if isfield(s_in.EXP, 'NMC'),
        NMC=s_in.EXP.NMC;
    elseif isfield(resolution, 'NMC'),
        NMC=resolution.NMC;
    else
        NMC  = 2000;
    end
else
    NMC = 2000;
    if isfield(s_in,'resolution')
        resolution = s_in.resolution;
    else
        error('No resolution field')
    end
end
%-----

% method: rescal5/rc_conv
% this code is very compact and efficient, after re-factoring and testing
% against rescal5.

if runGPU
    RMC = randn(4,NMC,'single','gpuArray');
else
    RMC = randn(4,NMC); % Monte Carlo points
end

% [rlu] [R] frame: Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV])
for frames={'rlu','spec'}  % others: 'cart','rlu_ABC','ABC'
    frame = resolution.(frames{1});
    M=frame.RM;
    [V,E]=eig(M);
    sigma=1./sqrt(diag(E)); % length along principal axes of gaussian
    
    % compute MC points on axes (with NMC points)
    xp   = bsxfun(@times,sigma,RMC);
    % get cloud in HKLE [rlu^3.meV], centred at 0.
    XMC  = inv(V)'*xp; % this is delta(HKLE) as a Gaussian distribution
    % compute the HKLE position in the lattice frame [a*,b*,c*,w]
    HKLE = frame.rlu2frame*resolution.HKLE(1:3)';
    HKLE(4)   = resolution.HKLE(4);
    HKLE = bsxfun(@plus,HKLE,XMC); % add HKLE location to Gaussian
    
    s_in.resolution.(frames{1}).cloud = HKLE;%{ HKLE(1,:)' HKLE(2,:)' HKLE(3,:)' HKLE(4,:)' }; % get 1D arrays per axis
    
end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% method: ResLib/ConvRes
% the code extracted from ConvRes does not seem to produce sensible Monte-Carlo
% points. The distribution is clearly not Gaussian, and extends very far from
% the HKLE position.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



% and now we can evaluate the function onto the axes.... and sum all values
% sum(feval(model, parameters, ax{:}))*resolution.R0/NMC

% the opposite operation (cloud -> RM) is computed from:
% <http://stackoverflow.com/questions/3417028/ellipse-around-the-data-in-matlab>
% B as columns
% Center = mean(B,1);
% X0 = bsxfun(@minus,B,Center);
% RM=X0'*X0 ./ (size(X0,1)-1);
