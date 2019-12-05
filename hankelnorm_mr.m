function [Ar, Br, Cr, Dr] = hankelnorm_mr(A,B,C,D,r)
%HANKELNORM_MR Creates a reduced model with respect to the hankel norm
    sys = ss(A,B,C,D);
    %% step 1: compute hsv and it's multiplicity
    hankel_sv = hsvd(sys);
    sv_r1 = hankel_sv(r+1);
    n = numel(hankel_sv);
    s = sum(abs(hankel_sv(r+1:end)-sv_r1)<(10*eps));
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
    U = B(r+1:r+s,:)*pinv(C(:,r+1:r+s)');
%     [U,~] = qr(U);
    %% step 5: Calculate all pass dilation
    Ab = inv(Gamma)*(sv_r1^2*A(1:end-s,1:end-s)'+P(1:end-s,1:end-s)*A(1:end-s,1:end-s)*P(1:end-s,1:end-s)-sv_r1*C(:,1:end-s)'*U*B(1:end-s,:)');
    Bb = inv(Gamma)*(P(1:end-s,1:end-s)*B(1:end-s,:)+sv_r1*C(:,1:end-s)'*U);
    Cb = C(:,1:end-s)*P(1:end-s,1:end-s)+sv_r1*U*B(1:end-s,:)';
    Db = D-sv_r1*U;
    %% step 6: sort based on eigenvalues
    [eigen_vec,eigen_val] = eig(Ab,'vector');
    is_stable = real(eigen_val)<0;
    n_stable = sum(is_stable);
    T_inv = [eigen_vec(:,is_stable) eigen_vec(:,~is_stable)];
    Ab = inv(T_inv)*Ab*T_inv;
    Bb = inv(T_inv)*Bb;
    Cb = Cb*T_inv;
    %% step 7: output reduced form
    Ar = Ab(1:n_stable,1:n_stable);
    Br = Bb(1:n_stable,:);
    Cr = Cb(:,1:n_stable);
    Dr = Db;
end
