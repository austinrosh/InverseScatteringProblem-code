%% detectorFields.m: 
% Computes the free-space scalar Green's Function for a given
% source/detector/wavenumber input.
function [ U_d ] = detectorFields(XX,YY,ZZ,xd,yd,zd,dim,k,Ndetect_x,Ndetect_y)

U_d = zeros(dim^3,(Ndetect_x)^2);

m = 1; %detector number
for i=1:Ndetect_x
    for j=1:Ndetect_y
        for p=1:(dim^3)
            r = sqrt(abs((XX(p)-xd(i))^2 + (YY(p)-yd(j))^2 + (ZZ(p)-zd)^2));
            U_d(p,m) = exp(-1j*k*r)/(4*pi*r);
        end
        m = m+1;
    end
    if m == Ndetect_x*Ndetect_y
        break
    end
end

end

