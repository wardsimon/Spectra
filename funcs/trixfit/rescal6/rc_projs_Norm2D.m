function [x1,y1,z1,x2,y2,z2]=rc_projs_Norm2D(R0,NP,i)
    %
    % MATLAB function to plot the projections of the resolution ellipse
    % of a triple axis
    %
    % DFM 10.11.95
        
    A=NP;
    % THIS IS FOR HWHM, but for our work, we are using sigma which is 1.
%     const=1.17741;
    const = 1;
    
    %----- Remove the vertical component from the matrix.
    
    B=[A(1,1:2),A(1,4);
        A(2,1:2),A(2,4);
        A(4,1:2),A(4,4)];
    
    %----- Work out projections for different cuts through the ellipse
    
    %----- S is the rotation matrix that diagonalises the projected ellipse
    
    if i ==1
            %----- 1. Qx, Qy plane
        [R0P,MP]=rc_int(3,R0,B);
        
        theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
        S=[];
        S=[cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        MP=S*MP*S';
        
        hwhm_xp=const/sqrt(MP(1,1));
        hwhm_yp=const/sqrt(MP(2,2));
        [x1,y1,z1]=rc_ellip2D(hwhm_xp,hwhm_yp,theta);
        
        %---------------- Add slice through Qx,Qy plane ----------------------
        
        MP=A(1:2,1:2);
        
        theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
        S=[];
        S=[cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        MP=S*MP*S';
        
        hwhm_xp=const/sqrt(MP(1,1));
        hwhm_yp=const/sqrt(MP(2,2));
        [x2,y2,z2]=rc_ellip2D(hwhm_xp,hwhm_yp,theta);
        
    elseif i==2
        %----- 2. Qx, W plane
        
        [R0P,MP]=rc_int(2,R0,B);
        
        theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
        S=[];
        S=[cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        MP=S*MP*S';
        
        hwhm_xp=const/sqrt(MP(1,1));
        hwhm_yp=const/sqrt(MP(2,2));        
        [x1,y1,z1]=rc_ellip2D(hwhm_xp,hwhm_yp,theta);

        
        %---------------- Add slice through Qx,W plane ----------------------
        
        MP=[A(1,1) A(1,4);A(4,1) A(4,4)];
        
        theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
        S=[];
        S=[cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        MP=S*MP*S';
        
        hwhm_xp=const/sqrt(MP(1,1));
        hwhm_yp=const/sqrt(MP(2,2));
        [x2,y2,z2]=rc_ellip2D(hwhm_xp,hwhm_yp,theta);
        
    elseif i==3
        %----- 3. Qy, W plane
        
        [R0P,MP]=rc_int(1,R0,B);
        
        theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
        S=[];
        S=[cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        MP=S*MP*S';
        
        hwhm_xp=const/sqrt(MP(1,1));
        hwhm_yp=const/sqrt(MP(2,2));
        [x1,y1,z1]=rc_ellip2D(hwhm_xp,hwhm_yp,theta);

        
        %---------------- Add slice through Qy,W plane ----------------------
        
        MP=[A(2,2) A(2,4);A(4,2) A(4,4)];
        
        theta=0.5*atan(2*MP(1,2)./(MP(1,1)-MP(2,2)));
        S=[];
        S=[cos(theta) sin(theta); -sin(theta) cos(theta)];
        
        MP=S*MP*S';
        
        hwhm_xp=const/sqrt(MP(1,1));
        hwhm_yp=const/sqrt(MP(2,2));
        [x2,y2,z2]=rc_ellip2D(hwhm_xp,hwhm_yp,theta);

    end
   
    if nargout == 5
        varargout{1} = theta;
    end
    
    
    
