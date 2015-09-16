function [p, dp, fixed,pnames]=mf_rpars
%
% MFIT function [p, dp, fixed,pnames]=mf_rpars
% 		Returns parameters read from parameters window
% 		MZ 29.11.94
%
[hmf_ctrl, hmf_data, hmf_pars]=mf_figs;

        p = [];
        dp = [];
        fixed = [];
        pnames = [];

if ~isempty(hmf_pars) & hmf_pars
        hm=get(hmf_pars,'userdata');

        [n,c]=size(hm);
          for i=1:n
                p = [ p str2num(get(hm(i,1),'String')) ];
                dp = [ dp str2num(get(hm(i,2),'String')) ];
                fixed = [ fixed get(hm(i,3),'Value') ];
                pnames = [ pnames ; get(hm(i,3),'string') ];
          end
end

if ~isempty(p) mf_upars(p,[]); end

