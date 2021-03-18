%% diag_Method
% this code implements the forced diagonal operator in order to compute the
% inverse recovery coefficients v1, v2, & v3 by solving for L1, L2, and L3
% in terms of forced-diagonal operations and known quantities

close
clear all
clc

format long

load forwardProblem.mat
load K1.mat
%load K2.mat

%load and declare forward operators
%K1 = K1_operator.K1;
%L1 = pinv(K1);
%K2 = K2_operator.K2;

%Load incident, detector, and scattered fields from the forward problem
U_i = FWD.U_i;
U_d = FWD.U_d;
U_s = FWD.U_s;
k = FWD.k;

%Load interaction matrix
G = FWD.G;

%Original material properties
v_original = FWD.V_vec;

%Notation used by Dr. Levinson
A = U_i;
B = U_d;
A_pi = pinv(U_i);
B_pi = pinv(U_d);

%% Compute recovery coefficients
v_hadd = pinv(((A'*A).*(B*B').')); %diag method operator

%first order
D_arg1 = diag(A'*U_s*B');
v1 = v_hadd*D_arg1;

%second order
v1_d = diag(v1);
D_arg2 = diag(A'*A*v1_d*G*v1_d*B*B');
v2 = -v_hadd*D_arg2;

%third order
v3_term1 = v_hadd*diag(A'*A*v1_d*G*v1_d*G*v1_d*B*B');
                                          
Q = v_hadd*diag(A'*A*v1_d*G*v1_d*B*B');   %Define P & Q as matricies contained in the other terms
Q_d = diag(Q);
P = v_hadd*diag(A'*A*v1_d*B*B');
P_d = diag(P);

v3_term2 = -v_hadd*diag(A'*A*P_d*G*Q_d*B*B');
v3_term3 = -v_hadd*diag(A'*A*Q_d*G*P_d*B*B');
v3 = -(v3_term1 + v3_term2 + v3_term3);



%% Recovery
v_firstOrder = v1;
v_secondOrder = v1 + v2; %should these be plus or minus
v_thirdOrder = v1 + v2 + v3;



error_firstOrder = norm(v_firstOrder-v_original)/norm(v_original)
error_secondOrder = norm(v_secondOrder-v_original)/norm(v_original)
error_thirdOrder = norm(v_thirdOrder-v_original)/norm(v_original)

dim = FWD.dim;
sample_csv1 = zeros(dim,dim);
sample_csv2 = zeros(dim,dim);
sample_csv3 = zeros(dim,dim);
sample_csv4 = zeros(dim,dim);
count = 1;
for i = 1:dim           %slice index
    while (count< i*dim^2)
        for j = 1:dim      %x i ndex
            for h = 1:dim   %y index
                sample_csv1(h,j) = abs(v_firstOrder(count));
                sample_csv2(h,j) = abs(v_secondOrder(count));
                sample_csv3(h,j) = abs(v_thirdOrder(count));
                original_cs(h,j) = abs(v_original(count));
                count = count+1;
            end
        end
    end
    
    figrows = ceil(dim/3);
    figcols = 3;
    figure(1)
    subplot(figrows,figcols,i)
    imagesc(sample_csv1)
    sgtitle('First Order Recovery')
    caxis([0 abs(max(v_firstOrder))])
    hold on
    title({'Slice: ', [num2str(i)]})
    %shading interp
    set(gca,'xticklabel',{[]}, 'yticklabel',{[]})
    colorbar
    
    figure(2)
    subplot(figrows,figcols,i)
    imagesc(sample_csv2)
    sgtitle('Second Order Recovery')
    caxis([0 abs(max(v_secondOrder))])
    hold on
    title({'Slice: ', [num2str(i)]})
   % shading interp
    set(gca,'xticklabel',{[]}, 'yticklabel',{[]})
    colorbar
    
    
    figure(3)
    subplot(figrows,figcols,i)
    imagesc(sample_csv3)
    sgtitle('Third Order Recovery')
    caxis([0 abs(max(v_thirdOrder))])
    hold on
    title({'Slice: ', [num2str(i)]})
    %shading interp
    set(gca,'xticklabel',{[]}, 'yticklabel',{[]})
    colorbar
    
    figure(4)
    subplot(figrows,figcols,i)
    imagesc(original_cs)
    sgtitle('Original')
    caxis([0 abs(max(v_original))])
    hold on
    title({'Slice: ', [num2str(i)]})
    %shading interp
    set(gca,'xticklabel',{[]}, 'yticklabel',{[]})
    colorbar
    
    count = i*dim^2;
end