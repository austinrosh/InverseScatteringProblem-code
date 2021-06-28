%% diag_Method(filename,flag)
% this code implements the forced diagonal operator in order to compute the
% inverse recovery coefficients v1, v2, & v3 by solving for L1, L2, and L3
% in terms of forced-diagonal operations and known quantities
%
% <filename> specifies the .mat file to load
% <flag> specifies whether or not to plot the recovery cross-sections. 0 =
% no plotting, 1 = plotting.

function [v1,v2,v3,v4] = diagMethod(filename,flag)
    format long
    load([filename,'.mat']) 

    %Load incident, detector, and scattered fields from the forward problem
    A = FWD.U_i;
    B = FWD.U_d;
    U_s = FWD.U_s;

    %Load scattering interaction matrix
    G = FWD.G;

    %Original material properties
    v_original = FWD.V_vec;
    
    %% Born Series Recovery Terms
    v_hadd = pinv(((A'*A).*(B*B').')); %diag method operator
    
    %% Compute recovery coefficients
    v_hadd = pinv(((A'*A).*(B*B').')); %diag method operator

    %first order (Born Approximation)
    D_arg1 = diag(A'*U_s*B');
    v1 = v_hadd*D_arg1; 
    v1_d = diag(v1);

    %% Born Series Recovery Terms
    
    %forward series terms
    %fixed point estimates from born series
    nu1 = A*v1_d*B;
    nu2 = A*v1_d*G*v1_d*B;
    nu3 = A*v1_d*G*v1_d*G*v1_d*B;
    nu4 = A*v1_d*G*v1_d*G*v1_d*G*v1_d*B;
    
    %inverse born series terms
    phi1 = v_hadd*diag(A'*(nu1)*B');
    phi1 = diag(phi1);
    phi2 = v_hadd*diag(A'*(nu2)*B');  
    phi2 = diag(phi2);
    
    phi3 = v_hadd*diag(A'*(nu3)*B');  
    phi3 = diag(phi3);
    phi4 = v_hadd*diag(A'*(nu4)*B');
    phi4 = diag(phi4);
   
    %% 
    %second order
    v1_d = diag(v1);
    D_arg2 = diag(A'*(A*v1_d*G*v1_d*B)*B');
    %D_arg2 = phi2;
    v2 = -v_hadd*D_arg2;
    %v2_d = diag(v2);
    v2_d = diag(v2);
    

    %third order
    v3_term1 = v_hadd*diag(A'*(A*v1_d*G*v1_d*G*v1_d*B)*B');
    %v3_term1 = phi3;
    v3_term2 = -v_hadd*diag(A'*(A*phi1*G*phi2*B)*B');
    v3_term3 = -v_hadd*diag(A'*(A*phi2*G*phi1*B)*B');
    v3 = -(v3_term1 + v3_term2 + v3_term3);

    %4th order
    %terms start to get mixed here
    v4_term1 = v_hadd*diag(A'*(A*v1_d*G*v1_d*G*v1_d*G*v1_d*B)*B');
    %v4_term1 = phi4;
    v4_term2 = v_hadd*diag(A'*(A*phi2*G*phi2*B)*B');
    
    v4_term3a = -v_hadd*diag(A'*(A*phi1*G*-diag((v_hadd*diag(A'*(A*v1_d*G*v2_d*B)*B')))*B)*B');
    v4_term3b = -v_hadd*diag(A'*(A*phi2*G*-diag((v_hadd*diag(A'*(A*v2_d*B)*B')))*B)*B');
    %v4_term3c = v_hadd*diag(A'*(A*v1_d*G*v1_d*G*-(v_hadd*diag(A'*(A*v2_d*B)*B')))*B');
    v4_term3c = v_hadd*diag(A'*(A*v1_d*G*v1_d*G*-v2_d*B)*B');
    v4_term3 = -(v4_term3a + v4_term3b + v4_term3c);
    
    v4_term4a =  -v_hadd*diag(A'*(A*phi1*G*-diag((v_hadd*diag(A'*(A*v2_d*G*v1_d*B)*B')))*B)*B');
    v4_term4b = -v_hadd*diag(A'*(A*-diag((v_hadd*diag(A'*(A*v1_d*G*v2_d*B)*B')))*G*phi1*B)*B');
    %v4_term4c = v_hadd*diag(A'*(A*v1_d*G*v1_d*G*-(v_hadd*diag(A'*(A*v2_d*B)*B')))*B');
    v4_term4c = v_hadd*diag(A'*(A*v1_d*G*-v2_d*G*v1_d*B)*B');
    v4_term4 = -(v4_term4a + v4_term4b + v4_term4c);
    
    v4_term5a =  -v_hadd*diag(A'*(A*-diag((v_hadd*diag(A'*(A*v2_d*B)*B')))*G*phi2*B)*B');
    v4_term5b = -v_hadd*diag(A'*(A*-diag((v_hadd*diag(A'*(A*v2_d*G*v1_d*B)*B')))*G*phi1*B)*B');
    %v4_term5c = v_hadd*diag(A'*(A*v1_d*G*v1_d*G*-(v_hadd*diag(A'*(A*v2_d*B)*B')))*B');
    v4_term5c = v_hadd*diag(A'*(A*-v2_d*G*v1_d*G*v1_d*B)*B');
    v4_term5 = -(v4_term5a + v4_term5b + v4_term5c);
    
    v4 = v4_term1 + v4_term2 + v4_term3 + v4_term4 + v4_term5;
    
   
    
    

    %% Recovery
    v_firstOrder = v1;
    v_secondOrder = v1 + v2; 
    v_thirdOrder = v1 + v2 + v3;
    v_fourthOrder = v1 + v2 + v3 + v4;

    error_firstOrder = norm(v_firstOrder-v_original)/norm(v_original);
    err1 = error_firstOrder;
    error_secondOrder = norm(v_secondOrder-v_original)/norm(v_original);
    err2 = error_secondOrder;
    error_thirdOrder = norm(v_thirdOrder-v_original)/norm(v_original);
    err3 = error_thirdOrder;
    error_fourthOrder = norm(v_fourthOrder-v_original)/norm(v_original);
    err4 = error_fourthOrder;
    
    display(err1)
    display(err2)
    display(err3)
    display(err4)

    %% Plotting

    dim = FWD.dim;
    sample_csv1 = zeros(dim,dim);
    sample_csv2 = zeros(dim,dim);
    sample_csv3 = zeros(dim,dim);
    sample_csv4 = zeros(dim,dim);
    count = 1;
    if(flag==1)
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

end
