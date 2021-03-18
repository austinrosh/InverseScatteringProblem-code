%% K2 Coefficient
% This code computes the K2 coefficient of the Born-Series expansion in
% order to solve for the material vector V using the inverse Born-Series
% expansion
clc
close
clear all

load forwardProblem.mat

x = FWD.V_vec;
b = FWD.U_s(:);

U_i = FWD.U_i;
U_d = FWD.U_d;
U_s = FWD.U_s;
G = FWD.G;

xd = FWD.xd;
yd = FWD.yd;
zd = FWD.zd;

x1 = FWD.V_vec;
x2 = x1;
%x = kron(x1,x2);
n = 1;
for i = 1:length(x1)
    for j = 1:length(x1)
        x(n) = x1(i)*x1(j);
        n = n+1;
    end
end

numSources = FWD.numSources;
numDetectors = FWD.numDetectors;
dim = FWD.dim;

%% Create K2 matrix
K2 = ones(numSources*numDetectors, (dim^3)^2);

m = 1; %source index
n = 1; %detector index
p = 1; %first medium index
h = 1; %second meidun index

for i = 1:numSources*numDetectors
    for j = 1:(dim^3)^2
        K2(i,j) = U_i(m,p)*G(p,h)*U_d(h,n);
        p = p+1; %go to next medium index
        if (p == dim^3 + 1) %if the entire first medium as been searched
            p = 1; %set first medium index back to 1
            h = h+1; %increment second medium index
        end
        if (h == dim^3 + 1) %if entire second medium has been searched
            h = 1; %set second medium index back to 1
        end
    end
    m = m+1; %increment source
    if (m == numSources + 1) %if we have considered all sources for the current detector
        m = 1; %set source index back to 1
        n = n+1; %increment detector index
    end
   
end
K2_operator.K2 = K2;
save('K2.mat','K2_operator')

%% Compare K2 solution w/ Ui*V*G*v*Ud
% Check that the solution via the K2 operator is the same as given by the
% Born series for the second order scattering term
U_s2coeff = K2*x;
V = diag(FWD.V_vec);
U_s2series = U_i*V*G*V*U_d;
U_s2series = U_s2series(:);
error = norm(U_s2coeff-U_s2series)/norm(U_s2series);