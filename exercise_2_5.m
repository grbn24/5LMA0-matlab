A = [0 1 0; 0 0 1; -1 2 3];
B = [0; 0; 1];
C = [0 1 2];
x = [1; 1; 1;];
t_max = 10;
P = zeros(3,3,t_max);
Q = zeros(3,3,t_max);
e_reach = zeros(t_max,1);
e_obs = zeros(t_max,1);
P(:,:,1) = B*B';
Q(:,:,1) = C'*C;
e_reach(1) = x'*pinv(P(:,:,1))*x;
e_obs(1) = x'*Q(:,:,1)*x;

for t = 2:t_max
    P(:,:,t) = A*P(:,:,t-1)*A'+B*B';
    Q(:,:,t) = A*Q(:,:,t-1)*A'+C'*C;
    e_reach(t) = x'*inv(P(:,:,t))*x;
    e_obs(t) = x'*Q(:,:,t)*x;
end
e_obs(e_obs>100) = NaN;
figure(1);
plot(3:t_max,e_reach(3:end));
hold on
% plot(3:t_max,e_obs);
hold off