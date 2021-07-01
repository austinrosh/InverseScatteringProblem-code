%% diag_Method.m
% this code implements the forced diagonal operator in order to compute the
% inverse recovery coefficients v_i by solving for the inverse series terms
% L_i in terms of forced-diagonal operations and known quantities
%
% <filename> specifies the .mat file to load
% <flag> specifies whether or not to plot the recovery cross-sections. 0 = no plotting, 1 = plotting.

clc
clear
close all
format long

filename = 'forwardProblem';
load([filename,'.mat']) 
flag = 0; %change to 1 to show plots in diagMethod() function

A = FWD.U_i; %incident fields
B = FWD.U_d; %detector fields
U_s = FWD.U_s; %measured scattered field
G = FWD.G; %scattering interaction matrix
v_hadd = pinv(((A'*A).*(B*B').')); %diag method operator
v_original = FWD.V_vec; %original material

%% V1 & Static Terms
%1st order (Born Approximation)
v1 = L1(U_s,FWD);
v1_d = diag(v1);

%declare static terms reused in subsequent recovery coefficients up to
%the Nth term to be computed

arg1 = A*v1_d*B; %note that this should be equal to U_s
arg2 = A*v1_d*G*v1_d*B;
arg3 = A*v1_d*G*v1_d*G*v1_d*B;
arg4 = A*v1_d*G*v1_d*G*v1_d*G*v1_d*B;
arg5 = A*v1_d*G*v1_d*G*v1_d*G*v1_d*G*v1_d*B;
arg6 = A*v1_d*G*v1_d*G*v1_d*G*v1_d*G*v1_d*G*v1_d*B;


%% ISP Terms 2...N
% N = 6 in this code   
%2nd order
v2 = L2(arg1,arg1,FWD);

%3rd order
v3 = L3(arg1,arg2,arg3,FWD);
%v3_d = diag(v3);
%4th order

v4 = L4(arg1,arg2,arg3,arg4,FWD);

 %5th order
v5 = L5(arg1,arg2,arg3,arg4,arg5,FWD);

% %6th order
% v6 = L6(nu1,nu2,nu3,nu4,nu5,arg6,FWD);

v_ISP = [v1 v2 v3 v4 v5 v6];
%% Recovery
v_firstOrder = v1;
v_secondOrder = v1 + v2; 
v_thirdOrder = v1 + v2 + v3;
v_fourthOrder = v1 + v2 + v3 + v4;
v_fifthOrder = v1 + v2 + v3 + v4 + v5;
v_sixthOrder = v1 + v2 + v3 + v4 + v5 + v6;

%error differences
error_firstOrder = norm(v_firstOrder-v_original)/norm(v_original);
err1 = error_firstOrder;
error_secondOrder = norm(v_secondOrder-v_original)/norm(v_original);
err2 = error_secondOrder;
error_thirdOrder = norm(v_thirdOrder-v_original)/norm(v_original);
err3 = error_thirdOrder;
error_fourthOrder = norm(v_fourthOrder-v_original)/norm(v_original);
err4 = error_fourthOrder;
error_fifthOrder = norm(v_fifthOrder-v_original)/norm(v_original);
err5 = error_fifthOrder;
error_sixthOrder = norm(v_sixthOrder-v_original)/norm(v_original);
err6 = error_sixthOrder;

display(err1)
display(err2)
display(err3)
display(err4)
display(err5)
display(err6)


%% Plotting
%{
if(flag==1)   
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
                    sample_csv4(h,j) = abs(v_fourthOrder(count));
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
        %caxis([0 abs(max(v_secondOrder))])
        caxis([0 abs(max(v_firstOrder))])

        hold on
        title({'Slice: ', [num2str(i)]})
       % shading interp
        set(gca,'xticklabel',{[]}, 'yticklabel',{[]})
        colorbar


        figure(3)
        subplot(figrows,figcols,i)
        imagesc(sample_csv3)
        sgtitle('Third Order Recovery')
        %caxis([0 abs(max(v_thirdOrder))])
        caxis([0 abs(max(v_firstOrder))])

        hold on
        title({'Slice: ', [num2str(i)]})
        %shading interp
        set(gca,'xticklabel',{[]}, 'yticklabel',{[]})
        colorbar

        figure(4)
        subplot(figrows,figcols,i)
        imagesc(sample_csv4)
        sgtitle('Fourth Order Recovery')
        %caxis([0 abs(max(v_fourthOrder))])
        caxis([0 abs(max(v_firstOrder))])

        hold on
        title({'Slice: ', [num2str(i)]})
        %shading interp
        set(gca,'xticklabel',{[]}, 'yticklabel',{[]})
        colorbar


        figure(5)
        subplot(figrows,figcols,i)
        imagesc(original_cs)
        sgtitle('Original')
        %caxis([0 abs(max(v_original))])
        caxis([0 abs(max(v_firstOrder))])

        hold on
        title({'Slice: ', [num2str(i)]})
        %shading interp
        set(gca,'xticklabel',{[]}, 'yticklabel',{[]})
        colorbar


        count = i*dim^2;
    end
end
%}


%% Recovery Term Functions

%Xi is a forward scattering coeffieicnt given by (A)*(material1)*(G)*(material2)*G...*(materialN)*(B)
%Xi = Phi_i

function v = L1(X1,FWD)
    A = FWD.U_i; %incident fields
    B = FWD.U_d; %detector fields
    v_hadd = pinv(((A'*A).*(B*B').')); %diag method operator
    v = v_hadd*diag(A'*X1*B'); 
end

function v = L2(X1,X2,FWD)
    A = FWD.U_i; %incident fields
    B = FWD.U_d; %detector fields
    G = FWD.G; %scattering interaction matrix
    M1 = diag(L1(X1,FWD));
    M2 = diag(L1(X2,FWD));
    Y1 = A*M1*G*M2*B; %scattering term derived from X1 & X2
    v = -L1(Y1,FWD);
end

function v = L3(X1,X2,X3,FWD)
    Y1 = L1(X3,FWD);
    Y2 = L2(X1,X2,FWD)+L2(X2,X1,FWD);
    v = -(Y1+Y2);
end

function v = L4(X1,X2,X3,X4,FWD)
    Y1 = L1(X4,FWD);
    Y2 = L2(X2,X2,FWD) + L2(X3,X1,FWD) + L2(X1,X3,FWD);
    Y3 = L3(X1,X1,X2,FWD) + L3(X1,X2,X1,FWD) + L3(X2,X1,X1,FWD);
    v = -(Y1+Y2+Y3);
end

function v = L5(X1,X2,X3,X4,X5,FWD)
    Y1 = L1(X5,FWD);
    Y2 = L2(X1,X4,FWD)+L2(X4,X1,FWD)+L2(X2,X3,FWD)+L2(X3,X2,FWD);
    Y3 = L3(X1,X1,X3,FWD)+L3(X1,X3,X1,FWD)+L3(X3,X1,X1,FWD)+L3(X1,X2,X2,FWD)+L3(X2,X1,X2,FWD)+L3(X2,X2,X1,FWD);
    Y4 = L4(X1,X1,X1,X2,FWD)+L4(X1,X1,X2,X1,FWD)+L4(X1,X2,X1,X1,FWD)+L4(X2,X1,X1,X1,FWD);
    v = -(Y1+Y2+Y3+Y4);
end

function v = L6(X1,X2,X3,X4,X5,X6,FWD)
    Y1 = L1(X6,FWD);
    Y2 = L2(X3,X3,FWD)+L2(X1,X5,FWD)+L2(X5,X1,FWD)+L2(X2,X4,FWD)+L2(X4,X2,FWD);
    Y3 = L3(X2,X2,X2,FWD)+L3(X1,X1,X4,FWD)+L3(X1,X4,X1,FWD)+L3(X4,X1,X1,FWD)+...
        L3(X1,X2,X3,FWD)+L3(X1,X3,X2,FWD)+L3(X2,X1,X3,FWD)+L3(X2,X3,X1,FWD)+L3(X3,X1,X2,FWD)+L3(X3,X2,X1,FWD);
    Y4 = L4(X1,X1,X1,X3,FWD)+L4(X1,X1,X3,X1,FWD)+L4(X1,X3,X1,X1,FWD)+L4(X3,X1,X1,X1,FWD)+...
        L4(X1,X1,X2,X2,FWD)+L4(X1,X2,X1,X2,FWD)+L4(X1,X2,X2,X1,FWD)+L4(X2,X1,X1,X2,FWD)+L4(X2,X1,X2,X1,FWD)+L4(X2,X2,X1,X1,FWD);
    Y5 = L5(X1,X1,X1,X1,X2,FWD)+L5(X1,X1,X1,X2,X1,FWD)+L5(X1,X1,X2,X1,X1,FWD)+L5(X1,X2,X1,X1,X1,FWD)+L5(X2,X1,X1,X1,X1,FWD);
    v = -(Y1+Y2+Y3+Y4+Y5);
end
    
