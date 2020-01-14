Q = [2 0; 0 1];
W = [ 1 0 -1; 0 1 0];
Q_c = chol(Q);
[U, S, V] = svd(Q_c*W);
phi = inv(Q_c)*U;
