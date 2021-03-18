%% K1 Coefficient
% Computes the K1 coefficient in the inverse-born series expansion using
% the inverse operator neumann-series method. K1 is the Born Approximation
% recovery of the Born-Series.
clc
close
clear all

load forwardProblem.mat

x = FWD.V_vec;
b = FWD.U_s(:);

U_i = FWD.U_i;
U_d = FWD.U_d;
U_s = FWD.U_s;

xd = FWD.xd;
yd = FWD.yd;
zd = FWD.zd;

x = FWD.V_vec;
b = U_s(:);

numSources = FWD.numSources;
numDetectors = FWD.numDetectors;
dim = FWD.dim;
%% create K1 matrix
K1 = zeros(numSources*numDetectors,dim^3);
m = 1; %source index
n = 1; %detector index
for i = 1:numSources*numDetectors
    for j = 1:dim^3
        K1(i,j) = U_i(m,j)*U_d(j,n);
    end
    m = m+1;
    if (m == numSources+1)
        m = 1;
        n = n+1;
    end
    
%     if(n == numDetectors)
%         break
%     end
end

U_s1coeff = K1*x;

V = diag(FWD.V_vec);
U_s1series = U_i*V*U_d;
U_s1series = U_s1series(:);
error = norm(U_s1coeff-U_s1series)/norm(U_s1series)
K1_operator.K1 = K1;
save('K1.mat', 'K1_operator')

%other ways of recovering V from the scattered field measurements
%{
[u1,s1,v1] = svd(A);
s1_inv = pinv(s1); 
x_ls1 = v1*s1_inv*u1'*b; %compute LS solution via SVD
x_ls2 = pinv(A)*b; %compute LS solution by pseudo-inverse
[u2,s2,v2] = svd(U_i);
[u3,s3,v3] = svd(U_d);
x_ls3 = A\b;
%}

%---------------------------------------------------------------------%
%--------------------- Regularization & Imaging ----------------------%
%{
L = pinv(U_i);
R = pinv(U_d);
V = diag(U_s);
V_ls = L*U_s*R;
x_ls = diag(V_ls);
err1 = norm(x_ls-x);
disp('----------------------------------------------------')
fprintf('Error b/w recovered and true V values using direct pseudo-inverse = %f \n', err1);
%x_ = A\b;

%% Regularization
% Tikhonov regularization is applied in order to get a better condition
% number, which will allow for more accurate recovery of the material V by
% means of taking the pseudoinverse of A

lambda = 1e-6;
A_tik = (A'*A + (lambda^2).*eye(dim^3)); %Tikhonov regularization of the normal matrix used for SVD
x_tik = (A_tik)\(A'*b);
kappa = cond(A_tik);
err2 = norm(x-x_tik);
fprintf('Tikhonov regularized (lambda = %e) \n', lambda)
fprintf('Normal matrix condition number: %e \n', kappa);
fprintf('Error b/w recovered and true V using Tikhonov regularization: %f \n', err2)

%% Image the recovered solution
sample_cs = zeros(dim,dim);
%p = [1.*(1:dim)];
count = 1;
for i = 1:dim           %slice index
    while (count< i*dim^2)
        for j = 1:dim      %x i ndex
            for h = 1:dim   %y index
                sample_cs(h,j) = abs(x_ls(count));
                count = count+1;
            end
        end
    end
    subplot(5,3,i)
    imagesc(sample_cs)
    caxis([0 abs(max(x_ls))])
    hold on
    title({'Slice: ', [num2str(i)]})
    %shading interp
%     xticks([1 5 10 15])
%     xticklabels({'0','0.25','0.75','1'})
%     yticks([1 5 10 15])
%     yticklabels({'0','0.25','0.75','1'})
    set(gca,'xticklabel',{[]}, 'yticklabel',{[]})
    count = i*dim^2;
end
%}
    
        

