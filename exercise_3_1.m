load building.mat
A = full(A);
ev = eig(A);
logicalstr = ["false", "true"];
isstable = all(real(ev)<0);
reach = ctrb(A,B);
obs = obsv(A,C);
sys = ss(A,B,C,D);
[H,wout] = freqresp(sys);
%% Compute Gramians
P = gram(sys,'c');
Q = gram(sys,'o');
isreach = rank(P)==rank(A);
isobs = rank(Q)==rank(A);
fprintf("System stable: %s\nSystem reachable: %s\nSystem observable: %s\n",logicalstr(isstable+1),logicalstr(isreach+1),logicalstr(isobs+1));
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
r = 8;
Ar = An(1:r,1:r);
Br = Bn(1:r,:);
Cr = Cn(:,1:r);
Dr = D-Cn*inv(An)*Bn+Cr*inv(Ar)*Br;
sys_r = ss(Ar,Br,Cr,Dr);
Pr = gram(sys_r,'c');
Qr = gram(sys_r,'o');
[Hr,woutr] = freqresp(sys_r,wout);
isstable_r = all(real(eig(Ar))<0);
isreach_r = rank(Pr)==rank(Ar);
isobs_r = rank(Qr)==rank(Ar);
fprintf("Reduced system stable: %s\nReduced system reachable: %s\nReduced system observable: %s\n",logicalstr(isstable_r+1),logicalstr(isreach_r+1),logicalstr(isobs_r+1));

%% comparison plots
figure(1);
subplot(2,2,1);
loglog(wout,abs(H(:)))
hold on;
loglog(woutr,abs(Hr(:)));
title('Frequency Response');
xlabel('Frequency [rad/s]');
ylabel('Powah');
set(gca, 'XScale','log');
legend(["Normal", sprintf("Reduced r = %d",r)]);
hold off
subplot(2,2,2)
[y,t] = step(sys);
plot(t,y);
hold on
[y,t]=step(sys_r);
plot(t,y)
title("Step Response");
legend(["Normal", sprintf("Reduced r = %d",r)]);
hold off
%% upper bound on error
G = tf(sys);
G_r = tf(sys_r);
max_e = norm(G-G_r,inf);
fprintf("Maximum error %.5f\n",max_e);
subplot(2,2,3:4)
bode(G,G_r);
legend(["Normal", sprintf("Reduced r = %d",r)]);