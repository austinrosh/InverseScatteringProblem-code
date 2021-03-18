%% Forward Scattering Problem
% This code solves the scalar Helmholtz equation on a three-dimensional
% domain in the presence of a scattering object which interacts with an
% incident field taking the form of a plane wave. Multiple sources and
% detectors are able to be considered.

clc
close all
clear 
%% Initialize grid
Nx = 10;
hx = 1/Nx;
h = hx;

Ny = Nx;
hy = 1/Ny;
Nz = Nx;
hz = 1/Nz;

x = 0:1/Nx:1;
y = 0:1/Ny:1;
z = 0:1/Nz:1;

dim = length(x);

[X,Y,Z] = meshgrid(x,y,z);
XX = X(:);
YY = Y(:);
ZZ = Z(:);

samp_start = 0.2;
samp_end = 0.35;

rho = 4e-2;                    %permittivity of scatterer

lambda = 0.3;                     %wavelength
k = (2*pi)/lambda;                %wavenumber
%% Initialize sources (plane wave)
Ntheta = 5;                     %total of Ntheta*Nphi sources
Nphi = 5; 

U_i = createSources(XX,YY,ZZ,k,dim,Ntheta,Nphi);

%% Initialize material matrix
% V is a diagonal matrix with a defined material of k^2*eta in the sample
% region. The following lines of code initialize the material matrix to
% have vector elements which correspond to the desired X,Y, and Z coordinate
% combinations.
VV_vec = zeros(dim^3,1);

for j=1:(dim^3)
    if ( ((XX(j)>=samp_start)&&(YY(j)>=samp_start)&&(ZZ(j)>=samp_start))&&...
            ((XX(j)<=samp_end)&&(YY(j)<=samp_end)&&(ZZ(j)<=samp_end)))
   
                                 %only allow permissible region to have defined permittivity
          VV_vec(j) = (k^2)*rho; %norm(GV) determines convergence of Born series, must be less than 1 to ensure convergence
    end
end
VV = diag(VV_vec);          %VV is a diagonal matrix with elements corresponding to where the mateiral is defined in the grid space
%% Initialize detector locations and calculate fields at detector

%goal is to compute the scattered field at the detector location(s)

Ndetect_x= 5;              %numebr of x detectors
Ndetect_y= 5;              %number of y detectors
%Ndetectors = Ndetect_x*Ndetect_y; %total number of detectors
[xd,yd,zd] = createDetectors(Ndetect_x,Ndetect_y);
U_d = detectorFields(XX,YY,ZZ,xd,yd,zd,dim,k,Ndetect_x,Ndetect_y);

%% Compute Green's Function
%compute the Green's function relating the detector to the scatterer for
%each point in the sample, for each detector.
G = greensMatrix(XX,YY,ZZ,dim,k,h);
%G = (h^3).*G;
%% Compute scattered field

%compute the scattered field for every X & Y combination in the sample at
%the detector locations
I = eye(dim^3,dim^3);
U_s = U_i*((I-VV*G)\VV)*U_d;
%GV_norm = norm(VV*G,inf);

 %% Save forward problem data
VG = norm(VV*G);
fprintf('The norm |diag(V)*G| = %f \n', VG)
FWD.V_vec = VV_vec;
FWD.U_i = U_i;
FWD.U_s = U_s;
FWD.U_d = U_d;
FWD.G = G;
FWD.xd = xd;
FWD.yd = yd;
FWD.zd = zd;
FWD.k = k;
%FWD.fieldDim = Ndetect_x*Ndetect_y*Ntheta*Nphi;
FWD.dim = dim;
FWD.numSources = Ntheta*Nphi;
FWD.numDetectors = Ndetect_x*Ndetect_y;
save('forwardProblem.mat','FWD');