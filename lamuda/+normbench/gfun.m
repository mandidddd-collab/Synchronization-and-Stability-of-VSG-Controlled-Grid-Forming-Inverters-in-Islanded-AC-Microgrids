function value = gfun(x)
%GFUN g(x)=sqrt(x)-(x-1)/2 acos((x-1)/(x+1)), x>=1.
validateattributes(x,{'numeric'},{'real','finite','>=',1});
value = sqrt(x)-(x-1).*atan(1./sqrt(x));
end
