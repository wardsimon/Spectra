function varargout = rc_res()
    %
    % MATLAB function to calculate the resolution function
    % of a triple axis
    %
    % DFM 10.11.95
    
    
    %----- f converts from energy units into k^2
    %      f=0.4826 for meV and f=1.996854 for THz
    
    f=0.4826;
    
    %----- Get method from Parameter window
    
    method=get(findobj('tag','hrc_rescal_method'),'Userdata');
    
    %----- Get parameters from Parameter and Instrumentation window
    
    p=rc_savp('respars');
    
    %----- Calculate Q vector
    
    [Q2c,Qmag]=rc_re2rc([p(19) p(20) p(21)], [p(22) p(23) p(24)], [p(31) p(32) p(33)]);
    
    %----- Calculate resolution matrix in Q frame, and check that scattering triangle closes
    
    [R0,NP,vi,vf,Error]=feval(method,f,Qmag,p,0);
    
    if Error ~= 0; disp('Scattering triangle does not close'), return; end
    
    %----- Test of Normalisation constant
    
    [R0P,NPP]=rc_int(1,R0,NP);
    [R0P,NPP]=rc_int(1,R0P,NPP);
    [R0P,NPP]=rc_int(1,R0P,NPP);
    
    disp(' Resolution Matrix in frame Qx, Qy, Qz')
    
    NP
    
    %----- Calculate Bragg widths
    
    disp(' Bragg widths')
    
    [bragw]=rc_bragg(NP)
    
    %----- Diagonalise resolution matrix in Q frame
    
    [a,b]=eig(NP);
    
    a'*NP*a;
    
    inv(a)'*b*inv(a);
    
    %----- Now work out transformations
    
    Q2c;
    
    A1=[p(25) p(26) p(27)]';
    A2=[p(28) p(29) p(30)]';
    
    V1=Q2c*A1;
    V2=Q2c*A2;
    
    %----- Form unit vectors V1, V2, V3 in scattering plane
    
    V3=cross(V1,V2);
    V2=cross(V3,V1);
    V3=V3/sqrt(sum(V3.*V3));
    V2=V2/sqrt(sum(V2.*V2));
    V1=V1/sqrt(sum(V1.*V1));
    
    U=[V1';V2';V3'];
    
    %----- S transformation matrix from (h,k,l) to V1,V2,V3
    
    S=U*Q2c;
    S_save=S;
    
    %----- Work out angle of Q wrt to V1, V2
    
    TT=S*[p(31) p(32) p(33)]';
    cos_theta=TT(1)/sqrt(sum(TT.*TT));
    sin_theta=TT(2)/sqrt(sum(TT.*TT));
    
    %----- Rotation matrix from Q to V1,V2,V3
    
    R=[cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1];
    
    T=zeros(4,4);
    T(4,4)=1;
    T(1:3,1:3)=R*S;
    
    disp(' Resolution matrix in frame V1, V2, V3')
    
    V_norm = T'*NP*T
    
    %----- Phonon Widths
    
    %----- Transform scan and dispersion vectors into A-1 and meV vectors
    %      in the V1,V2,V3 coordinate frame.
    
    DQ_cart=(T*[p(35) p(36) p(37) p(38)]')';
    
    scan=DQ_cart/sqrt(sum(DQ_cart.*DQ_cart));
    plane=(T(1:3,1:3)*[p(39) p(40) p(41)]')';
    plane=[p(42)*plane/sqrt(sum(plane.*plane)) -1];
    plane=plane/sqrt(sum(plane.*plane));
    
    % ----- Calculate the phonon width.
    
    [phoni,phonw]=rc_phon(1,NP,plane);
    phonw=abs((phonw/(sum(scan.*plane))))
    
    %----- Diagonalised matrix in RLU
    
    %disp(' Resolution matrix in RLU and Energy')
    [V,E]=eig(T'*NP*T);
    
    %----- Projection axes in terms of reciprocal space axes ---------------
    
    qx=[p(31) p(32) p(33)]';
    qy=cross(inv(Q2c)*V3,qx);
    qx=qx/sqrt(sum(qx.*qx))
    qy=qy/sqrt(sum(qy.*qy))
    
    
    %----- Write parameters to figure
    
    rc_lab(bragw,phonw,R0,vi,vf,qx,qy);
    
    %----- Projections of the Resolution ellipse in Qx, Qy, Qz, W space
    
    rc_projs(R0,NP);
    
j =1;
        width(j,1) = 1/sqrt(E(1,1));
        width(j,2) = 1/sqrt(E(2,2));
        width(j,3) = 1/sqrt(E(3,3));
        width(j,4) = 1/sqrt(E(4,4));
    
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
    title(BigAx,sprintf('[ %0.03f %0.03f %0.03f %0.03f ] in r.l.u',p(31:34)))
    
    for m = 1:numel(AX)
        set(AX(m),'ButtonDownFcn',@hitme,'Userdata',sprintf('%i',m))
    end
    
    function hitme(varargin)
        [c, b] = ind2sub([4,4],str2double(get(varargin{1},'UserData')));
        fprintf('FWHM: %s\t= %0.05f\n',t{b},fwhm(b))
        fprintf('FWHM: %s\t= %0.05f\n',t{c},fwhm(c))
    end
end