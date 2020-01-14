function [sys_red] = match_max_analytical(sys, freq, order)
%MATCH_MAX_ANALYTICAL Summary of this function goes here
%   Detailed explanation goes here
    [A,B,C,D] = ssdata(sys);
    I = eye(size(A));
    D0 = D + C*inv(I*freq-A)*B;
    Ap = inv(I*freq-A);
    Bp = Ap*B;
    Cp = C*Ap;
    sys_p = ss(Ap,Bp,Cp,D0);
    N0 = gram(sys_p,'o');
    N0 = N0(1:order);
    R0 = gram(sys_p,'c');
    R0 = R0(1:order);
    H = N0*R0;
    R = chol(H);
    N = R';
    U = R0*inv(R);
    V = [inv(N)*N0]';
    A_r = V'*A*U;
    B_r = -V'*B;
    C_r = C*U;
    D_r = D;
    sys_red = ss(A_r,B_r,C_r,D_r);
end

