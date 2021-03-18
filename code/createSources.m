%% createSources.m
% creates a specified number of sources given by Ntheta and Nphi on the
% interval theta in [0,pi] and phi in [0,pi]
function [uu_i] = createSources(X,Y,Z,k,dim,Ntheta,Nphi)

theta = linspace(0,pi,Ntheta+1);
theta = theta(1:Ntheta);
phi = linspace(0,pi,Nphi);

uu_i = zeros(Ntheta*Nphi,dim^3);    %initialize source matrix with proper dimension;
                                    %each row is a different source
                                    %combination with unique theta/phi
                                    %angles
m = 1; %source index
for i=1:length(theta)
    for j=1:length(phi)
        kx = k*sin(theta(i))*cos(phi(j));
        ky = k*sin(theta(i))*sin(phi(j));
        kz = k*cos(theta(i));
        u_i = exp(-1i*(kx*X + ky*Y +kz*Z));
        uu_i(m,:) = u_i;
        m = m+1;
    end
end
end

