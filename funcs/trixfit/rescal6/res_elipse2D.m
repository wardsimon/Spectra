function [X1, X2, X3, Y1, Y2, Y3, Z1, Z2, Z3] = res_elipse2D(x,scan,parameter_file,configuration_file,varargin)
    global verbose
    if isempty(verbose)
        verbose = 0;
    end
    fid=fopen(configuration_file);
    if fid == -1
        error('CFG file not found');
    end
    ppop=parameter_read(fid);
    fid=fopen(parameter_file);
    if fid == -1
        error('PAR file not found');
    end
    pcn=parameter_read(fid);
    p=[pcn; ppop];
    
    R_Do = 0;
    f = 0.4826;
    
    method = 'rc_popma';
    
    points = zeros(length(x),4);
    points(:,scan(4)) = x;
    a = 1:4; a(scan(4))=[];
    points(:,a) = scan(ones(length(x),1),1:3);
    
    width = zeros(length(x),4); %Qh, Qk, Ql, E
    R0_corrected = zeros(length(x),1);
    NP_hkl = zeros(length(x),4,4);
    for j = 1:length(x)
        [Q2c] = rc_re2rc([p(19) p(20) p(21)], ...
            [p(22) p(23) p(24)],...
            points(j,1:3));
        % Q2c: matrix to transform a vector Q(H,K,L) into cart. coord. w. resp. to a
        %      basis with a || x and b in the x-y-plane
        
        %----- Now work out transformations
        
        % A1,A2: the rec. lattice vectors in the scattering plane
        A1 = [p(25) p(26) p(27)]';
        A2 = [p(28) p(29) p(30)]';
        
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
        
        p(31) = points(j,1);
        p(32) = points(j,2);
        p(33) = points(j,3);
        p(34) = points(j,4);
        
        %----- Calculate focusing curvatures
        % Note method has to contain rc_popma!
        if ~isempty(strfind(method,'rc_popma'))
            
            % NOTE! the scan might have manual curvatures!
            if R_Do == 1
                rho = rc_focus(p);
                % NOTE Rho is in [1/m]
                rmh = rho(1);% mon. vert. foc. in 1/m
                rmv = rho(2);% mon. vert. foc. in 1/m
                rah = rho(3);% ana. horz. foc. in 1/m
                rav = rho(4);% ana. vert. foc. in 1/m
            else
                rav = p(69);
                rah = p(68);
                rmv = p(67);
                rmh = p(66);
            end
            
            % We migth not have focusing optics on the monochromator
            % or analyser. So if given values are ~0 we say they are flat
            p(66) = (p(66)>1E-5)*rmh;
            p(67) = (p(67)>1E-5)*rmv;
            p(68) = (p(68)>1E-5)*rah;
            p(69) = (p(69)>1E-5)*rav;
            % rc_popma converts these to cm^-1 i.e. pres/100
        end
        
        %----- Q vector in cartesian coordinates
        Qcart = Q2c*p(31:33);
        Qmag  = sqrt(sum(Qcart.*Qcart));
        
        %----- Resolution matrix in Q frame. NP is resolution matrix in Qx, Qy & Qz
        [R0,NP,vi,vf,Error] = feval(method,f,Qmag,p',1);
        
        %----- Use correction of resolution volume
        
        R0_corrected(j) = R0/sqrt(det(NP))/(2*pi)^2;
        
        %----- Work out angle of Q wrt to V1, V2
        TT = S*[p(31) p(32) p(33)]';% Q in V1,V2,V3 coord. sys.
        cos_theta = TT(1)/sqrt(sum(TT.*TT));
        sin_theta = TT(2)/sqrt(sum(TT.*TT));
        
        %----- Rotation matrix from Q to V1,V2,V3
        
        R = [cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1];
        T = zeros(4,4);
        T(4,4) = 1;
        T(1:3,1:3) = R*S;
        
        %----- Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV]) i.e. V1, V2, V3 frame!
        
        NP_hkl(j,:,:) = T'*NP*T; % This is the matrix of the quadratic form NP, when the Q-vectors are given
        % in (h,k,l) coordinates
        try
            [V,E] = eig(reshape(NP_hkl(j,:,:),4,4));% V contains the normalized eigenvectors of M,
        catch
            error('Q point [%0.03f %0.03f %0.03f %0.03f] could not be calculated',pres(31),pres(32),pres(33),pres(34))
        end
        width(j,1) = 1/sqrt(E(1,1));
        width(j,2) = 1/sqrt(E(2,2));
        width(j,3) = 1/sqrt(E(3,3));
        width(j,4) = 1/sqrt(E(4,4));
        %         b_mat = reshape(inv((V'))',1,16);
        
        
        %         fwhm = width(j,:);
        % create the Monte Carlo points with the correct
        % standard deviation
        %         xp = fwhm(ones(NMC,1),:)'.*randn(4,NMC);
        %
        %         XMC = reshape(b_mat(1:16),4,4)'*xp;
        %
        %         Qx = XMC(A1==1,:) + points(j,A1==1);            % 1 x NMC matrix
        %         Qy = XMC(A2==1,:) + points(j,A2==1);            % 1 x NMC matrix
        %         En = XMC(4,:)     + points(j,4);                % 1 x NMC matrix
        
        % We need to make NP_hkl into the correct form. A1, A2, A3 EN where A3 is out of plane!
        % This might not work for A = [1 1 0] etc!
        neps = 1E-10;
        A = squeeze(NP_hkl(j,:,:));
        A(:,[find(abs(A(1,:))<neps) 3]) = A(:,[3 find(abs(A(1,:))<neps)]);
        A([find(abs(A(:,1))<neps) 3],:) = A([3 find(abs(A(:,1))<neps)],:);
        
        % We have a normal distribution, lets say to 3 sigma i.e. 99.73% of all points!
        if isempty(varargin)
            cont_percent = 0.9973;
            n_sig = erfinv(cont_percent)*sqrt(2);
        else
            if isnumeric(varargin{1})
                n_sig = varargin{1};
            else
                cont_percent = 0.9973;
                n_sig = erfinv(cont_percent)*sqrt(2);
            end
        end
        
        [X1, Y1, Z1, ~, ~] = rc_projs_Norm2D(R0_corrected,A,1);
        [X2, Y2, Z2, ~, ~] = rc_projs_Norm2D(R0_corrected,A,2);
        [X3, Y3, Z3, ~, ~] = rc_projs_Norm2D(R0_corrected,A,3);
        
        X1 = X1*n_sig + points(j,abs(A1)==1);
        X2 = X2*n_sig + points(j,abs(A1)==1);
        X3 = X3*n_sig + points(j,abs(A2)==1);
        
        Y1 = Y1*n_sig + points(j,abs(A2)==1);
        Y2 = Y2*n_sig + points(j,4);
        Y3 = Y3*n_sig + points(j,4);
        
        
        if verbose
            NMC = 1E4;
            fwhm = width(j,:);
            xp = fwhm(ones(NMC,1),:)'.*randn(4,NMC);
            XMC = (V')\xp;
            
           figure
           [H,AX,BigAx,P,PAx] = plotmatrix(XMC');
           t{1} = sprintf('[%i %i %i]',A1(:));
           t{2} = sprintf('[%i %i %i]',A2(:));
           t{3} = sprintf('[%i %i %i]',cross(A1,A2));
           t{4} =         'Energy';
           for m = 1:length(t)
               ylabel(AX(m,1),t{m})
               xlabel(AX(4,m),t{m})
           end
           title(BigAx,sprintf('[ %0.03f %0.03f %0.03f %0.03f ]',points(j,:)))
           
           for m = 1:numel(AX)
               set(AX(m),'ButtonDownFcn',@hitme,'Userdata',sprintf('%i',m))
           end
           
        end
        
    end
    function hitme(varargin)
        [c, b] = ind2sub([4,4],str2double(get(varargin{1},'UserData')));
        fprintf('FWHM: %s\t= %0.05f\n',t{b},fwhm(b))
        fprintf('FWHM: %s\t= %0.05f\n',t{c},fwhm(c))
    end
end

%=======================================================
function [p]=parameter_read(fid)
    %
    %
    %
    %-------------- Initialize arrays-----------------
    
    data=[];
    header='';
    text=fgetl(fid);
    
    %----------------- Load data-----------------------
    
    while (text>0)
        [temp count]=sscanf(text,'%f');
        if isempty(temp)
            header=[header text];
        else
            if (count==size(data,2) | isempty(data))
                data=[data; temp'];
            end
        end
        text=fgetl(fid);
    end
    fclose(fid);
    
    p=data(:,1);
end