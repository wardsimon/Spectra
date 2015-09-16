function [y, name, pnames, pin]=user(x,p, flag)
% user      : user defined constant for Multifunc
% [y, {name, pnames, pin}]= user(x,p, {flag})
%
% MFIT constant
% parameter name is asked when flag = 'user_ask'
% to specify pname : flag = pname

% Author:  EF
% Description:  user parameter

if nargin==2;
	y=p(1)*ones(size(x));
else
	y=[];
	name='User Constant';
	pnames=str2mat('UserConstant');
	if flag==1, pin=1; else pin = p; end
	if flag==2
		mf_msg('Click on user constant(y axis)');
		[x1 y1]=ginput(1);
		pin=y1;
	end
	if ischar(flag)
		if strcmp(flag,'user_ask')
			prompt = inputdlg('Enter parameter name','User parameter');
		else
			prompt = flag;
		end
		if iscell(prompt)
			prompt = char(prompt{1});
		end
		if ~isempty(prompt)
			name = prompt;
			pnames = prompt;
		end
	end
end
