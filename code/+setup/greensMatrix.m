%% greensMatrix.m
% computes the interaction between all points in the sample region with a
% diven discretization, incident wave-number, and geometric dimension.
function [G] = greensMatrix(XX,YY,ZZ,dim,k,h)

G = zeros(dim^3,dim^3);
zeta = log(26+15*sqrt(3))-pi/2; %singularity correction

for i=1:dim^3
    for j=1:dim^3
        r = sqrt(abs(XX(j)-XX(i))^2 + (YY(j)-YY(i))^2 + (ZZ(j)-ZZ(i))^2);
        if (r == 0)
            G(i,j) = ((h)^2)*(zeta-1j*k*h); %singularity treatment
        else
            G(i,j) = (h^3)*(exp(-1j*k*r)./(4.*pi.*r));
        end
    end
end

end
