function [val] = gfunc(x,y,z,kval)
    r   = sqrt(x.*x+y.*y+z.*z);
    val = -exp(1i*kval*r)./(4*pi*r);
end

