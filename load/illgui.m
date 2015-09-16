function [x,y,err,xlab,ylab,monitor]=illgui(filename)
%
%function [x,y,err,xlab,ylab,monitor]=illgui(filename)
%
%   This function extracts scans from data files of the type generated at
%   the ILL using the new format. Search scan variable and ask user for 
%   the X,Y variables to import. EF 08.07.97 DFM 12.6.96
%
%Output variable
%   x        - the independent variable
%   y        - the dependent variable
%   err      - uncertainty in dependent variable = sqrt(y)
%   xlab     - name of 'x' data
%   ylab     - name of 'y' data

% uses : ffind.c as a mex file. illbatch.m

[x, y, err, xlab, ylab,monitor]=illbatch([ filename ',gui' ]);




