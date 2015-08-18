function [x,y,err,xlab,ylab,monitor]=tasgui(filename)
%
%function [x,y,err,xlab,ylab,monitor]=tasgui(filename)
%
%   This function extracts scans from TASCOM data files
%   (old or  new format). Ask user for 
%   the X,Y variables to import. EF 08.07.97 DFM 12.6.96
%
%Output variable
%   x        - the independent variable
%   y        - the dependent variable
%   err      - uncertainty in dependent variable = sqrt(y)
%   xlab     - name of 'x' data
%   ylab     - name of 'y' data

% uses : ffind.c as a mex file. illbatch.m

[x, y, err, xlab, ylab,monitor]=tasbatch([ filename ',gui' ]);




