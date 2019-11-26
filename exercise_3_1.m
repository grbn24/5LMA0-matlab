load building.mat
A = full(A);
ev = eig(A);
logicalstr = ["false", "true"];
isstable = all(real(ev)<0);
reach = ctrb(A,B);
obs = obsv(A,C);
isreach = rank(reach)==rank(A);
isobs = rank(obs)==rank(A);
fprintf("System stable: %s\nSystem reachable: %s\nSystem observable: %s\n",logicalstr(isstable+1),logicalstr(isreach+1),logicalstr(isobs+1));
sys = ss(A,B,C,D);
[H,wout] = freqresp(sys);
figure(1);
subplot(2,2,1);
loglog(wout,abs(H(:)))
title('Frequency Response');
xlabel('Frequency [rad/s]');
ylabel('Powah');
subplot(2,2,2);
step(sys);
%% Compute Gramians
P = gram(sys,'c');
Q = gram(sys,'o');
%% Make balanced representation
R = chol(P); 
R = R'; % matlab choleski suxx
[U,S,~] = svd(R'*Q*R);
T = (S.^0.25)*U'*inv(R);
An = T*A*inv(T);
Bn = T*B;
Cn = C*inv(T);
Dn = D;
sys_n = ss(An,Bn,Cn,Dn);
Pn = gram(sys_n,'c');
Qn = gram(sys_n,'o');
%% make reduced system
r = 5;
Ar = An(1:r,1:r);
Br = Bn(1:r,:);
Cr = Cn(:,1:r);
Dr = D-Cn*inv(An)*Bn+Cr*inv(Ar)*Br;
sys_r = ss(Ar,Br,Cr,Dr);
[Hr,woutr] = freqresp(sys_r);
%% comparison plots
subplot(2,2,3)
loglog(woutr,abs(Hr(:)));
title(sprintf('frequency response reduce r = %d',r));
xlabel('Frequency [rad/s]')
ylabel('Powah')
subplot(2,2,4)
step(sys_r)
title(sprintf('Step Response reduced r = %d',r));
%% upper bound on error
G = tf(sys);
G_r = tf(sys_r);
max_e = norm(G-G_r,inf);
fprintf("Maximum error %.5f\n",max_e);