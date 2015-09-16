function [X, Y, Z] = rc_ellip2D(sigma_x,sigma_y,theta)
    
    % Confidence region. We need this if we have artifacts.
    ext = 4*erfinv(0.999)*sqrt(2);
    
    % To normalise to volume or not
    A = 1/(2*pi*sigma_x*sigma_y);
    A = 1;
    % This needs a clockwise rotation. We have an anti-clockwise rotation 
    theta = -theta;
    x = linspace((-sigma_x*ext),(sigma_x*ext),500);
    y = linspace((-sigma_y*ext),(sigma_y*ext),500);
    a = cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
    b = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2 ;
    c = sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;
    
    [X, Y] = meshgrid(x,y);
    Z = A*exp( - (a*X.^2 + 2*b*X.*Y + c*Y.^2)) ;
end