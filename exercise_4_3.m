n_vec = [1 5 10 50];
Ts = -1;
for n = n_vec
    if n > 0
        A = [zeros(n-1,1) eye(n-1);zeros(1,n)];
        B = [zeros(n-1,1); 1];
        C = [1 zeros(1,n-1)];
        D = [0];
    elseif n == 0
        A = 0;
        B = 0;
        C = 0;
        D = 1;
    end
    sys = ss(A,B,C,D,Ts);
    Q = dlyap(A',C'*C);
    P = dlyap(A,B*B');
    singular_val = svd(P*Q);
    hankel_norm = max(singular_val);
    fprintf("Max hankel singular value for n = %d is %.4f\n",n,hankel_norm);
end
r = 10;
sys_r = hankelmr(sys,r);
impulse(sys,sys_r)
