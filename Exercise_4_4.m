load('building.mat');
sys = ss(A,B,C,D);
h_sv = hsvd(sys);
h_sv_norm = h_sv./sum(h_sv);
r = find(cumsum(h_sv_norm)>0.95,1,'first');
sys_r = hankelmr(sys,r+1);
[Ar, Br, Cr, ~] = ssdata(sys_r);
Dr = D - C*inv(A)*B + Cr*inv(Ar)*Br;
sys_r = ss(Ar,Br,Cr,Dr);
figure(1);
subplot(2,1,1)
step(sys,sys_r);
subplot(2,1,2)
bode(sys,sys_r);
h_norm = max(hsvd(tf(sys)-tf(sys_r)))