function s_out = lt(s1,s2)
%
% function r = lt(s1,s2)
%
% @SPEC1D/lt function to give the product of spectrum s1 with a value or spectrum s2.
%
% Simon Ward 26/01/2016 - simon.ward@psi.ch
%

if isa(s1,'spec1d')
    s = s1;
    val = s2;
else
    s = s2;
    val = s1;
end

for i = 1:length(s)
   ind = s(i).x < val;
   r = s(i);
   r.x = s(i).x(ind);
   r.y = s(i).y(ind);
   r.e = s(i).e(ind);
   if ~isempty(r.yfit)
       r.yfit = s(i).yfit(ind);
   end
   s_out(i) = feval(class(r),r);
end