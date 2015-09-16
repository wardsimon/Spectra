function [y, name, pnames, pin]=heaviside(x, p)
% heaviside : Heaviside function
% MFIT Heaviside fitting function
% p = [ Amp Centre FullWidth BackG ]


if nargin==2
    y = p(4) + (abs(x - p(2)) < p(3)/2)*(p(1) - p(4));
else
        y=[];
        name='heaviside';
        pnames=str2mat('Amplitude','Centre','Width','Background');
        if flag==1, pin=[1 0 1 1]; else pin = p; end
        if flag==2
                mf_msg('Click on peak');
                [cen amp]=ginput(1);
                mf_msg('Click on width');
                [width y]=ginput(1);
                width=abs(width-cen);
                mf_msg('Click on background');
                [x bg]=ginput(1);
                amp=amp-bg;
                pin=[amp cen width bg];
        end
end
