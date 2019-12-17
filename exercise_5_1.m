%% inputs
r = 1;
N = 10;

%% statestpace creation and simulation
A = [0.5 0; 1 -0.5];
B = [4;0];
C = [1 1];
D = 0;
sys_x = ss(A,B,C,D,1);
sys_state = sys_x;
sys_state.C = eye(2);
[X] = lsim(sys_state,[1 0 0 0], 0:3, [0 0])'; % X should have n_states rows and n_t columns
%% create POD basis
[U,S,~] = svd(X);
phi1 = U(:,1);





% BULLSHIT
% Ar = U(:,1:r)'*A*U(:,1:r);
% Br = U(:,1:r)'*B;
% Cr = eye(2)*U(:,1:r);
% Dr = D - eye(2)*inv(A)*B + Cr*inv(Ar)*Br;
% sys_r = ss(Ar,Br, Cr,Dr,1);
% subplot(1,3,1)
% step(sys_state)
% subplot(1,3,2)
% step(sys_r,'r')
% subplot(1,3,3)
% step(sys_state,sys_r)