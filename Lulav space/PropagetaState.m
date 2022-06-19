function [Lambda,Omega,t] = PropagetaState(Lambda,Omega,t,new_t)

invA = StateEq(t,new_t);



Lambda = invA'*Lambda*invA;
Omega = invA'*Omega;
t = new_t;