function mf_upars(p, dp, fixed)
%
% MFIT function mf_upars(p, dp, fixed)
% Update parameters displayed in parameter window
% MZ 29.11.94
%
if nargin < 2
  dp=[];
end
if nargin  < 3
  fixed = [];
end

upars = findobj('tag','mf_ParWindow');
if isempty(upars)
  mf_newfn;
  upars = findobj('tag','mf_ParWindow');
  disp('Open ParWin')
end

h=get(upars(1),'userdata');

if isempty(h)
    return
end

[npars,i] = size(h);
% npars=max ([ length(p) length(dp) length(fixed) ]);
if isempty(find(npars == [ length(p) length(dp) length(fixed) ]))
  disp('Parameter number is wrong')
  return;
end

for i=1:npars
  if ~isempty(p)
    set(h(i,1),'string',num2str(p(i),6));
  end
    if ~isempty(dp)
      set(h(i,2),'string',num2str(dp(i),6));
  end
    if ~isempty(fixed)
      set(h(i,3),'value',(fixed(i) ~= 0));
  end
end



