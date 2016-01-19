function [y, name, pnames, pin] = ResCalFit(x,p,flag)
    %% ResCalFit function [y, name, pnames, pin] = ResCalFit(x,p,flag)
    % For fitting triple-axis data using ResLibCal as a back-end.
    % S(Q,w) is convoluted with the 4D resolution function.
    %
    % NOTE!!! This only works with @spec1d/fits version > 4.2
    %
    % S. Ward (simon.ward@psi.ch), September 2015
    % Version 1.3 (18.01.2016) - Support for GPU
    
    global ResFitScn iter_l multifit_ind
    persistent iter_mem runGPU
    
    if isempty(iter_mem)
        iter_mem = 1;
    end
    if isempty(runGPU)
        runGPU = sdext.getpref('gpuArray').val;
    end
    
    if nargin == 2;
        % Check to see if the fit has been initialised
        if isempty(ResFitScn)
            error('The fitting needs to be intitialised with ResCalFit_ini');
        end
        % Check if we are doing a multifit.
        if isempty(multifit_ind)
            multifit_ind = 1:length(x);
        end
        % Check that the inititalised fit is the same as the fit
        if length(x) ~= length(ResFitScn(multifit_ind))
            error('The initialised fit is not the same as the fit. Check ResCalFit_ini');
        end
        
        % What frame are we in? Only recompute every scan, not point! (Errors will crop up otherwise...)
        if iter_l ~= iter_mem
            ResFitScn(multifit_ind) = arrayfun(@updateMC,ResFitScn(multifit_ind));
%             for i = multifit_ind
%                 ResFitScn(i).resolution = updateMC(ResFitScn(i).EXP,ResFitScn(i).resolution);
%             end
            iter_mem = iter_l;
        end
        
        % Start the fitting
        y = zeros(size(x));
        j = 1;
        % Each point
        for i = multifit_ind
            % Get the MC points and reshape them to observed Q
            %Q = [(ResFitScn(i).resolution.abc.hkl2Frame\[ResFitScn(i).resolution.abc.cloud{1:3}]')' ResFitScn(i).resolution.abc.cloud{4}];
            Q = ResFitScn(i).resolution.rlu.cloud';
            if iscell(Q)
                Q = ([ResFitScn(i).resolution.rlu.cloud{:}]');
            end
            if isfield(ResFitScn(i).EXP,'NMC')
                NMC = ResFitScn(i).EXP.NMC;
            else
                if iscell(ResFitScn(i).resolution.rlu.cloud)
                    NMC = length(ResFitScn(i).resolution.rlu.cloud{1});
                else
                    NMC = length(ResFitScn(i).resolution.rlu.cloud(1,:));
                end
                ResFitScn(i).EXP.NMC = NMC;
            end
            % Calculate the cross section
            S = ResFitScn(i).resolution.R0*sum(feval(ResFitScn(i).SQW,Q,p))/NMC;
            % Add on a background.
            if runGPU
                y(j) = gather(S) + feval(ResFitScn(i).BG,ResFitScn(i).resolution.HKLE,p);
            else
                y(j) = S + feval(ResFitScn(i).BG,ResFitScn(i).resolution.HKLE,p);
            end
            j = j+1;
        end
    else
        if flag == 2
            % Start the fitting
            % Calculate the cross section
            Q = zeros(length(ResFitScn),4);
            for i = 1:length(ResFitScn)
                Q(i,:) = ResFitScn(i).resolution.HKLE;
            end
            S = feval(ResFitScn(1).SQW,Q,p);
            % Add on a background.
            y = S + feval(ResFitScn(i).BG,Q,p);
        else
            
            %----- Parameter names
            y = [];
            pin = [];
            name = 'ResCalFit';
            pnames = cell(length(p),1);
            for i=1:length(p)
                pnames{i} = ['p' num2str(i)];
            end
        end
        
    end
    
end