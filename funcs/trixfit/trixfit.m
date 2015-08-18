function [y, name, pnames, pin] = trixfit(x,p,flag)
    %
    % TRIXFIT function [y, name, pnames, pin]=trixfit(x,p,Q1,Q2, flag)
    % for fitting triple-axis data.
    %
    % S(Q,w) is convoluted with the 4D resolution function.
    %
    % Des McMorrow, October 2001
    % S. Ward, July 2010 - 8.08.13, simon.ward@psi.ch
    
    %----- Define global variables
    global global_trixfit
    global first_call
    global global_b_mat global_fwhm global_R0 Qvec
    global load_the_matrices
    %----- Units: f energy units into k^2 (0.4826 for meV, 1.996854 for THz)
    %      At the moment, only works for meV!
    
    if nargin == 2
        
        %----- Unpack parameters from global_trixfit
        %        p
        
        f = 0.4826;
        
        %-----  mon_flag: 1=monitor, 0=time
        
        mon_flag = global_trixfit.monitor_flag;
        
        %----- Monte Carlo Steps for integration
        
        NMC = global_trixfit.monte_carlo_samples;
        
        %----- Cross-section
        
        xsec = global_trixfit.xsec_file;
        
        %----- Background definition
        
        bkgd_file = global_trixfit.bkgd_file;
        
        %----- Correction file
        
        corr_file = global_trixfit.corr_file;
        
        %-----  Method
        
        method = global_trixfit.resolution_method;
        
        %----- Initialise y
        if global_trixfit.oversample>0
            xx = linspace(x(1),x(end),global_trixfit.oversample*length(x));
            yy = zeros(size(xx));
        else
            xx = x;
            yy = zeros(size(xx));
        end
        %----- Rescal params
        
        pres = global_trixfit.pres(:);
        
        %----- information needed to reconstruct the scan points for a non-equidistant scan:
        
        scan_len = length(xx); % number of points in scan
        if p(10) == 0
            global Q
            xvec = Q(:);
        else
            xvec = zeros(4,1);
            xvec(p(10)) = 1; % xvec is a vector in (Q,w) space with a non-zero entry in dir. of x (x= x-value of s1d-obj)
        end
        Qs = reshape(p(1:4),4,1);
        Qs(xvec == 1) = xx(1);
        Qe = reshape(p(5:8),4,1);
        Qe(xvec == 1) = xx(end);
        dQ = Qe-Qs;
        
        % Try to set first_call to 1 when change in Qvec or scan length
        if isempty(Qvec) % The first time we are called
            first_call = 1;
        else
            if first_call ~= 1 % Have we deliberatly set
                if size(Qvec,1) == scan_len % Have scan lengths changed
                    if (sum(Qvec(1,:)' == Qs) == 4) && (sum(Qvec(end,:)'== Qe) == 4) % Have Q vectors changed
                        first_call = 0;
                    else
                        first_call = 1;
                    end
                end
            end
        end
        
        %----- Calculate Q2c matrix
        if first_call == 1
            %             persistent_b_mat=zeros(scan_len,16);
            %             persistent_fwhm=zeros(scan_len,4);
            %             persistent_R0=zeros(scan_len,1);
            if license('test', 'Distrib_Computing_Toolbox') && (global_trixfit.parallel==1)
                if matlabpool('size') > 0
                    er = lasterror; %#ok<LERR>
                    if ~isempty(strfind(er.identifier,'parfor'))
                        matlabpool('close')
                        matlabpool('open')
                    end
                    % If not assume everything is OK
                else
                    matlabpool('open')
                end
                c = onCleanup(@()matlabpool('close'));
            end
            %
            % pres(19:21): lattice parameters a,b,c
            % pres(22:24): angles alpha, beta, gamma
            % pres(31:33): Q-point in [r.l.u.] coordinates
            [Q2c] = rc_re2rc([pres(19) pres(20) pres(21)], ...
                             [pres(22) pres(23) pres(24)],...
                             [pres(31) pres(32) pres(33)]);
            % Q2c: matrix to transform a vector Q(H,K,L) into cart. coord. w. resp. to a
            %      basis with a || x and b in the x-y-plane
            
            %----- Now work out transformations
            
            % A1,A2: the rec. lattice vectors in the scattering plane
            A1 = [pres(25) pres(26) pres(27)]';
            A2 = [pres(28) pres(29) pres(30)]';
            
            V1 = Q2c*A1;
            V2 = Q2c*A2;
            
            %----- Form unit vectors V1, V2, V3 in scattering plane
            
            V3 = cross(V1,V2);
            V2 = cross(V3,V1);
            V3 = V3/sqrt(sum(V3.*V3));
            V2 = V2/sqrt(sum(V2.*V2));
            V1 = V1/sqrt(sum(V1.*V1));
            
            U = [V1';V2';V3']; % transforms from orthonormal (Qx,Qy,Qz) basis to an orthonormal basis
            % (V1,V2,V3) with a* �� V1, V2 perp. a*, c* in the
            % (V1,V2)-plane, U has been verfied to be orthogonal
            
            %----- S transformation matrix from (h,k,l) to V1,V2,V3
            
            S = U*Q2c;     % This is used to bring the CN matrix into a defined frame.
            
            %----- Calculate resolution widths for scan etc
            
            for j = 1:scan_len;
                
                % pres(31:34): H,K,L, w of scan, p(1:4): H,K,L, w as given in mfit
                t = (repmat(xx(j),4,1)-Qs.*xvec)./(dQ.*xvec);
                t(isinf(t)) = 0; t(isnan(t)) = 0;
                pres(31) = Qs(1)+t(1)*dQ(1);
                pres(32) = Qs(2)+t(2)*dQ(2);
                pres(33) = Qs(3)+t(3)*dQ(3);
                pres(34) = Qs(4)+t(4)*dQ(4);
                
                %----- Calculate focusing curvatures
                % Note method has to contain rc_popma!
                if ~isempty(strfind(method,'rc_popma'))
                    
                    % NOTE! the scan might have manual curvatures!
                    if global_trixfit.R_Do == 1
                        rho = rc_focus(pres);
                        % NOTE Rho is in [1/m]
                        rmh = rho(1);% mon. vert. foc. in 1/m
                        rmv = rho(2);% mon. vert. foc. in 1/m
                        rah = rho(3);% ana. horz. foc. in 1/m
                        rav = rho(4);% ana. vert. foc. in 1/m
                    else
                        rav = pres(69);
                        rah = pres(68);
                        rmv = pres(67);
                        rmh = pres(66);
                    end
                    
                    % We migth not have focusing optics on the monochromator
                    % or analyser. So if given values are ~0 we say they are flat
                    pres(66) = (pres(66)>1E-5)*rmh;
                    pres(67) = (pres(67)>1E-5)*rmv;
                    pres(68) = (pres(68)>1E-5)*rah;
                    pres(69) = (pres(69)>1E-5)*rav;
                    % rc_popma converts these to cm^-1 i.e. pres/100
                end
                
                %----- Q vector in cartesian coordinates
                Qcart = Q2c*pres(31:33);
                Qmag  = sqrt(sum(Qcart.*Qcart));
                
                %----- Resolution matrix in Q frame. NP is resolution matrix in Qx, Qy & Qz
                [R0,NP,vi,vf,Error] = feval(method,f,Qmag,pres',mon_flag);
                
                %----- Use correction of resolution volume
                
                %     R0_corrected=R0/(sqrt(det(NP))/((2*pi)^2)); %this is contained in
                %       the orig. version (rc_int), but seems wrong when comparing it to reslib
                %      R0_corrected=R0*(sqrt(det(NP))/((2*pi)^2));
                %          [R0 R0_corrected (sqrt(det(NP))/((2*pi)^2))]
                % when applying Bertrand`s R0, all factors are included in the rc_popma:
                R0_corrected = R0/sqrt(det(NP))/(2*pi)^2;
                
                %----- Work out angle of Q wrt to V1, V2
                TT = S*[pres(31) pres(32) pres(33)]';% Q in V1,V2,V3 coord. sys.
                cos_theta = TT(1)/sqrt(sum(TT.*TT));
                sin_theta = TT(2)/sqrt(sum(TT.*TT));
                
                %----- Rotation matrix from Q to V1,V2,V3
                
                R = [cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1];
                T = zeros(4,4);
                T(4,4) = 1;
                T(1:3,1:3) = R*S;
                
                %----- Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV])
                
                M = T'*NP*T; % This is the matrix of the quadratic form NP, when the Q-vectors are given
                % in (h,k,l) coordinates
                try
                    [V,E] = eig(M);% V contains the normalized eigenvectors of M,
                catch
                    error('Q point [%0.03f %0.03f %0.03f %0.03f] could not be calculated',pres(31),pres(32),pres(33),pres(34))
                end
                
                fwhm = 1./sqrt(diag(E));
                
                b_mat = reshape(inv((V'))',1,16);
                
                global_b_mat(j,:) = b_mat;
                global_fwhm(j,:) = fwhm;
                global_R0(j) = R0_corrected;
            end
            first_call = 0;
            
        end
        
        %----- Calculate function
        if isempty(load_the_matrices)
            lm = 1;
        else
            lm = load_the_matrices;
        end
        
        if global_trixfit.parallel
            p_temp = p(11:(end-nargin(bkgd_file)));
            parfor j = 1:scan_len;
                %----- Evaluate the cross-section and do the Monte Carlo integration
                t = (repmat(xx(j),4,1)-Qs.*xvec)./(dQ.*xvec);
                t(isinf(t)) = 0;
                t(isnan(t)) = 0;
                temp = [Qs(1)+t(1)*dQ(1) Qs(2)+t(2)*dQ(2) Qs(3)+t(3)*dQ(3) Qs(4)+t(4)*dQ(4)];
                % We might have a mex file!
                if exist(xsec,'file') == 3
                    s = feval(xsec,NMC,temp,p_temp,global_fwhm(j,1:4),global_b_mat(j,1:16));
                else
                    s = feval(xsec,NMC,temp,p,global_fwhm(j,1:4),global_b_mat(j,1:16),lm);
                end                %----- Calculate normalised intensity
                yy(j) = global_R0(j)*s;
            end
            for j = 1:scan_len;
                t = (xx(j)-Qs'*xvec)/(dQ'*xvec);
                Qvec(j,:) = [Qs(1)+t*dQ(1) Qs(2)+t*dQ(2) Qs(3)+t*dQ(3) Qs(4)+t*dQ(4)];
            end
        else
            for j = 1:scan_len;
                %----- Evaluate the cross-section and do the Monte Carlo integration
                t = (repmat(xx(j),4,1)-Qs.*xvec)./(dQ.*xvec);
                t(isinf(t)) = 0;
                t(isnan(t)) = 0;
                temp = [Qs(1)+t(1)*dQ(1) Qs(2)+t(2)*dQ(2) Qs(3)+t(3)*dQ(3) Qs(4)+t(4)*dQ(4)];
                Qvec(j,:) = temp;
                % Check for a compiled cross section
                if exist(xsec,'file') == 3
                    s = feval(xsec,NMC,temp,p(11:(end-nargin(bkgd_file))),global_fwhm(j,1:4),global_b_mat(j,1:16));
                else
                    s = feval(xsec,NMC,temp,p,global_fwhm(j,1:4),global_b_mat(j,1:16),lm);
                end
                %----- Calculate normalised intensity
                yy(j) = global_R0(j)*s;
            end
        end
        
        %----- Apply correction factor to calculated intensity, eg for lambda/2 in
        %      incident beam
        yy = feval(corr_file,xx,yy);
        
        %----- Add background
        yy = yy+feval(bkgd_file,xx,p);
        
        try
            if global_trixfit.low_pass_filter >= 1
                if global_trixfit.low_pass_filter == 1
                    s_over = ceil(length(yy)/25);
                else
                    s_over = round(global_trixfit.low_pass_filter);
                end
                if license('test', 'Signal_Toolbox')
                    yy = medfilt1(yy,s_over);
                    % This produces an ugly shift if s_over is even!
                    if ~mod(s_over,2)
                        yy = interp1(xx-mean(diff(xx))/2,yy,xx,'pchip');
                    end
                else
                    yy = smooth(real(yy),s_over,'sgolay');
                end
            end
        catch
            if global_trixfit.low_pass_filter >= 1
                warning('trixfit:smoothing','Smoothing has failed, possibly not enough points.')
            end
        end
        
        if global_trixfit.oversample > 0
            y = interp1(xx,yy,x);
        else
            y = yy;
        end
        
        
    else
        
        %----- Parameter names
        y = [];
        pin = [];
        pnam_file = global_trixfit.pnam_file;
        name = 'trixfit';
        pnames = feval(pnam_file,p);
        
    end