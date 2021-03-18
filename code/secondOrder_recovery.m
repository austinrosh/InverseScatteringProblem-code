%% secondOrder_recovery
% This code computes the second order solution of the inverse Born series
% for recovering the scattering potential term V. The operators computed in
% K1_Coefficient.m and K2_Coefficient.m are used here to compute the
% inverse series. The errors associated with this material recovery and the
% original material values are studied.
clc
close
%clear all
format long

load forwardProblem.mat
load K1.mat
load K2.mat

%load and declare forward operators
K1 = K1_operator.K1;
K2 = K2_operator.K2;

Ui = FWD.U_i;
Ud = FWD.U_d;
Us = FWD.U_s;

dim = FWD.dim;
%% Compute solution

%first and second inverse terms
L1 = pinv(K1);
L1_tp = kron(L1,L1);
L2 = -L1*K2*L1_tp;
Us_tp = kron(Us,Us);

v1 = L1*Us(:);
v1_tp = kron(v1,v1);
x_2ndOrder = v1 - L1*K2*v1_tp(:);
x_original = FWD.V_vec;

secondOrder_error = norm(x_2ndOrder-x_original)%/norm(x_original)
firstOrder_error = norm(v1-x_original)%/norm(x_original)

%% Regularization
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
%}
%% Image the recovered solution
sample_csOG = zeros(dim,dim);
sample_csREC = zeros(dim,dim);
count = 1;
for i = 1:dim           %slice index
    while (count< i*dim^2)
        for j = 1:dim      %x i ndex
            for h = 1:dim   %y index
                sample_csOG(h,j) = abs(x_original(count));
                sample_csREC(h,j) = abs(x_2ndOrder(count));
                count = count+1;
            end
        end
    end
    figure(1)
    subplot(3,2,i)
    imagesc(sample_csOG)
    title('Original')
    caxis([0 abs(max(x_original))])
    hold on
    title({'Slice: ', [num2str(i)]})
    shading interp
    set(gca,'xticklabel',{[]}, 'yticklabel',{[]})
    
    figure(2)
    subplot(3,2,i)
    imagesc(sample_csREC)
    title('Recovered')
    caxis([0 abs(max(x_2ndOrder))])
    hold on
    title({'Slice: ', [num2str(i)]})
    shading interp
    set(gca,'xticklabel',{[]}, 'yticklabel',{[]})
    
    count = i*dim^2;
end
%}
