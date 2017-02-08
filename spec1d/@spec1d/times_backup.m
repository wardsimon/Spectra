function rvec=times(s1vec,s2vec)
%
% function r=times(s1,s2)
%
% @SPEC1D/TIMES Times for 1D spectra
%
% Valid expressions: 
%
% 1. r=s1.*s2      , where s1 and s2 are spectra
%							length(s1)==1 or length(s2)==1 or length(s1)==length(s2)
% 2. r=yc.*s1      , where yc multiplies the y axis
%							length(y)==1 or length(y)==length(s1)
% 3. r=[xc,yc].*s1 , where xc (yc) multiplies the x (y) axis
%							as above, [xc,yc]' is allowed.
%
% DFM 1.4.98, HMR 18.4.2001
%

if isa(s1vec,'spec1d') & isa(s2vec,'spec1d')
  if length(s1vec)==length(s2vec)
    n1=1:length(s1vec);  
    n2=n1;
  elseif length(s1vec)>length(s2vec)
    n1=1:length(s1vec);  
    n2=ones(size(n1));
  else
    n2=1:length(s2vec);  
    n1=ones(size(n2));
  end
  for nvec=1:length(n1)
    s1=s1vec(n1(nvec));
    s2=s2vec(n2(nvec));
    if length(s1.x)~=length(s2.x)
      disp('Spectra are not the same length')
      disp('Uses interpolation for second spectrum')
      s2=interpolate(s2,s1);
    end
    s1.e=sqrt((s1.y.*s2.e).^2+(s2.y.*s1.e).^2);
    s1.y=s1.y.*s2.y;
    s1.yfit=[];
    rvec(nvec)=s1;
  end
  return
elseif ~isa(s1vec,'spec1d')
   stmp=s1vec;
   s1vec=s2vec;
   s2vec=stmp;
end   

if ~isa(s2vec,'spec1d')
  if size(s2vec,1)==length(s1vec)
    fac=s2vec;
  elseif size(s2vec,2)==length(s1vec)
    fac=s2vec';
  elseif length(s2vec)==1
    fac=s2vec*ones(length(s1vec),1);
  else
    fac=[s2vec(1)*ones(length(s1vec),1) s2vec(2)*ones(length(s1vec),1)];
  end
  if size(fac,2)==1
    fac=[1+0*fac fac];
  end
  for nvec=1:length(s1vec)
    s1=s1vec(nvec);
    if size(fac,2)~=2
      error('Something went wrong - must debug')
    end
    s1.x=s1.x*fac(nvec,1);
    s1.e=s1.e*fac(nvec,2);
    s1.y=s1.y*fac(nvec,2);
    s1.yfit=s1.yfit*fac(nvec,2);
    rvec(nvec)=s1;
  end
else
   error('Something went wrong - must debug')
end