function [y, name, pnames, pin] = ResCalFit(x,p,flag)
    %% ResCalFit function [y, name, pnames, pin] = ResCalFit(x,p,flag)
    % For fitting triple-axis data using ResLibCal as a back-end.
    % S(Q,w) is convoluted with the 4D resolution function.
    %
    % NOTE!!! This only works with @spec1d/fits version > 4.2
    %
    % S. Ward (simon.ward@psi.ch), September 2015
    % Version 1.2 (14.09.2015)
    
    global ResFitScn iter_l multifit_ind
    persistent iter_mem
    
    if isempty(iter_mem)
        iter_mem = 1;
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
            for i=multifit_ind
                ResFitScn(i) = ResLibCal(ResFitScn(i),'script');
            end
            iter_mem = iter_l;
        end
        
        % Start the fitting
        y = zeros(size(x));
        j = 1;
        % Each point
        for i = multifit_ind
            % Get the MC points and reshape them to observed Q
            %Q = [(ResFitScn(i).resolution.abc.hkl2Frame\[ResFitScn(i).resolution.abc.cloud{1:3}]')' ResFitScn(i).resolution.abc.cloud{4}];
            Q = cell2mat(ResFitScn(i).resolution.rlu.cloud);
            NMC = length(ResFitScn(i).resolution.rlu.cloud{1});
            % Calculate the cross section
            S = sum(feval(ResFitScn(i).SQW,Q,p))/NMC;
            % Add on a background.
            y(j) = S + feval(ResFitScn(i).BG,ResFitScn(i).resolution.HKLE,p);
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