function [y, name, pnames, pin]=expon(x,p, flag)
% expon     : Exponential decay
% [y, {name, pnames, pin}]=expon(x,p, {flag}) 
%
% MFIT Exponential fitting function
% p = [ Amp tau bg ]

% Author:  EF <manuf@ldv.univ-montp2.fr>, MZ <mzinkin@sghms.ac.uk>
% Description: Exponential decay

if nargin==2;
    y=p(3)+p(1)*exp(-x/p(2));
else
        y=[];
        name='Exponential';
        pnames=str2mat('Amplitude','Tau','Background');
        if flag==1, pin=[0 1 1]; else pin = p; end
        if flag==2
                mf_msg('Click on background');
                [cen bg]=ginput(1);
                mf_msg('Click on point 1');
                [x1 y1]=ginput(1);
                mf_msg('Click on point 2');
                [x2 y2]=ginput(1);
      y1=y1-bg; y2=y2-bg;
      len=(x2-x1)/log(y1/y2);
      amp=y1/exp(-x1/len);
                pin=[amp len bg];
        end
end
