function [Lambda,Omega] = Fuse(Lambda,Omega,A,H,W,R,z)


Lambda = Lambda  + A'*H'*W'*inv(R)*W*H*A;
Omega = Omega  + A'*H'*W'*inv(R)*W*z;