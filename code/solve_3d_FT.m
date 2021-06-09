clc
clear
close all

kval = 5.1;
h = 2*pi/kval/8;
L = 7;

hcorr = 3.900264920001956;

xs = 0:h:L;
xs = [-xs(end:-1:2),xs];
ys = xs;
zs = xs;
n = numel(xs);

[X,Y,Z] = meshgrid(xs,ys,zs);


Lmax = ceil(log((2*L/h)+1)/log(2)); %what is Lmax?
Lmax = 2^Lmax;
xs_lrg = h*(0:Lmax);
xs_lrg = [xs_lrg,-xs_lrg(end-1:-1:2)]; 
ys_lrg = xs_lrg;
zs_lrg = xs_lrg;

[Xx,Yy,Zz] = meshgrid(xs_lrg,ys_lrg,zs_lrg);
Gf = setup.gfunc(Xx,Yy,Zz,kval);
Gf = Gf*h^3;
Gf(1,1,1) = -hcorr*h^2/(4*pi); %is this for singularity? why only on the first Green's function matrix of the entire array?
Gf = fftn(Gf);

vpot = exp(-((X-4).*(X-4)+(Y).*(Y)+(Z-3).*(Z-3))*6);
%vpot = sin(X).*cos(Y);

% now apply
vft = fftn(vpot,size(Gf)); %pads vpot to have trailing zeroes to allign dimension to Gf
vft = vft.*Gf; %Fourier convolution theorem: FT of the convolution of two functions is the pointwise product of the FT of the two functions.
vft = ifftn(vft);
vapplied = vft(1:n,1:n,1:n);

v2plot = squeeze(vapplied(round(n/2),:,:));
%v2plot = squeeze(vapplied(17,:,:));
pcolor(imag(v2plot))
