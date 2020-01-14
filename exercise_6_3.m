load('building.mat');
r = 10;
f = 10;
sys_r = match_max_analytical(ss(A,B,C,D),f,r);
bode(ss(A,B,C,D),sys_r)