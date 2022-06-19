function [x,P,t] = StateEstimate(Lambda,Omega,t_state,t)
numerical_epsilon = 1e-3;

if(min(eig(Lambda))<numerical_epsilon)
    P = Inf(size(Lambda));
    x = zeros(size(Omega));
    return;
end

Tr = diag(diag(Lambda).^-0.5);
P = Tr*inv(Tr*Lambda*Tr)*Tr;
% if(min(diag(P))<numerical_epsilon)
%     P = Inf(size(Lambda));
%     x = zeros(size(Omega));
%     return;
% end

x = P*Omega;

A = StateEq(t,t_state);
P = A*P*A';
x = A*x;

end