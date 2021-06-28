function [v] = nTerms_ISP(filename, numTerms)
    assert(numTerms>=1, '<numTerms> must pass a non-empty array')
    load([filename,'.mat'])
    A = FWD.U_i;
    B = FWD.U_d;
    G = FWD.G;
    U_s = FWD.U_s;
    dim = FWD.dim;
    v_hadd = pinv(((A'*A).*(B*B').')); %diag method operator
    v = zeros(dim^3, numTerms);
    v(:,1) = L1(U_s,A,B,v_hadd);
    v(:,2) = L2(v(:,1),v(:,1),A,B,G,v_hadd);
    %v(:,3) = L3(v(:,1),v(:,1),v(:,1),A,B,G,v_hadd);    
end

function [v1_] = L1(y1,A,B,v_hadd)
    v1_ = v_hadd*diag(A'*y1*B'); 
end

function [v2_] = L2(y1,y2,A,B,G,v_hadd)
    P = diag(L1(A*diag(y1)*B,A,B,v_hadd));
    Q = diag(L1(A*diag(y2)*B,A,B,v_hadd));
    v2_ = -v_hadd*diag((A'*A)*P*G*Q*(B*B'));
end

function [v3_] = L3(y1,y2,y3,A,B,G,v_hadd)
    P = v_hadd*diag((A'*A)*diag(y1)*G*diag(y1)*(B*B')); %first term
    Q = diag(v_hadd*diag((A'*A)*diag(y1)*G*diag(y2)*(B*B')));
    W = diag(v_hadd*diag((A'*A)*diag(y2)*G*diag(y1)*(B*B')));
    v3_ = -(P + Q + W);
end
    