load('building.mat');
sys = ss(A,B,C,D);
h_sv = hsvd(sys);
h_sv_norm = h_sv./sum(h_sv);
r = find(cumsum(h_sv_norm)>0.95,1,'first');
% r = 10;
sys_r = hankelmr(sys,r+1);
[Ar2, Br2, Cr2, Dr2,] = hankelnorm_mr(A,B,C,D,r);
sys_r2 = ss(Ar2,Br2,Cr2,Dr2);
[Ar, Br, Cr, ~] = ssdata(sys_r);
Dr = D - C*inv(A)*B + Cr*inv(Ar)*Br;
sys_r = ss(Ar,Br,Cr,Dr);
figure(1);
subplot(3,1,1)
step(sys,sys_r,sys_r2);
subplot(3,1,[2 3])
bode(sys,sys_r,sys_r2);
figure(2)
subplot(1,3,1)
hsvd(sys)
subplot(1,3,2)
hsvd(sys_r)
subplot(1,3,3)
hsvd(sys_r2)
h_norm = max(hsvd(tf(sys)-tf(sys_r)));