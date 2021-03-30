function [v] = nTerms_ISP(filename, numTerms)

assert(numTerms>=1, '<numTerms> must be an integer than or equal to 1')
 
load([filename,'.mat'])
A = FWD.U_i;
B = FWD.U_d;
G = FWD.G;
U_s = FWD.U_s;
v_hadd = pinv(((A'*A).*(B*B').')); %diag method operator

%declare inverse terms
%need to update method to declare empty struct for speed & clarity
v = [];

%recursively find all inverse terms
N = numTerms;
v = recursiveDiag(v,A,B,G,U_s,v_hadd,N)

end

function [v] = recursiveDiag(v,A,B,G,U_s,v_hadd,N)

    if N == 1
        v{1} = v_hadd*diag(A'*U_s*B');
    else
        v = pi;
    end
end