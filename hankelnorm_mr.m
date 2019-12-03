function [Ar, Br, Cr, Dr] = hankelnorm_mr(A,B,C,D,r)
%HANKELNORM_MR Creates a reduced model with respect to the hankel norm
    sys = statespace(A,B,C,D);
    %% step 1: compute hsv and it's multiplicity
    hankel_sv = hsvd(sys);
    sv_r1 = hankel_sv(r+1);
    n = numel(hankel_sv);
    s = count(abs(hankel_sv(r+1:end)-sv_r1)<(10*eps));
    %% step 2: make structured P and Q
    P = zeros(n);
    P(1:n-s,1:n-s) = diag([hankel_sv(1:r);hankel_sv(r+s+1:end)]);
    P(n-s+1:end,n-s+1:end) = eye(s)*sv_r1;
    Q = P;
    %% step 3: make system in the same structure and define gamma
    Ap = A([1:r r+s+1:end r+1:r+s],[1:r r+s+1:end r+1:r+s]);
    Bp = B([1:r r+s+1:end r+1:r+s],:);
    Cp = C(:,[1:r r+s+1:end r+1:r+s]);
    Dp = D;
    Gamma = P(1:n-s,1:n-s)^2-sv_r1^2*eye(n-s);
    %% step 4: Get unitary U
    
end

