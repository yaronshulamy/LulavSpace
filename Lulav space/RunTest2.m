

alpha0 =10;%m/sec^2
%noise params
rou_b_a = 1e-3;
rou_k_a = 1e-5;
rou_k_p = 1e-1;
t1= 100;%sec
p0=0;%m
v0=0;%m/sec
freq1 = 100;%Hz
delay = 0.5; %sec
freq2 = 1;%Hz

pdd= -alpha0;%m/sec^2

t=0:(1/freq1):t1;

p = p0 + t.*v0+ t.^2/2*pdd;
v =  v0 + t.*pdd;






b_a = randn(1)*rou_b_a;
k_a = randn(size(t))*rou_k_a;


am = pdd + b_a + k_a;

p_mask = (freq1/freq2/2):(freq1/freq2):length(t); 

tpm = t(p_mask);

k_p = randn(size(tpm))*rou_k_p;

pm = p(p_mask) + k_p;

%x = [p,v,a,b];
t_state = t(end);
E = 0;
Omega = zeros(4,1);
Lambda = zeros(4,4);
jjj=1;
figure
for iii=1:length(t)
    dt = t(iii) - t_state ;
    A = [1,dt,dt^2/2,0;0,1,dt,0;0,0,1,0;0,0,0,1];
    H = [0,0,1,1];
    R = rou_k_a^2;
    W =1;
    %E = E + (H*A*x-am(iii))'*W'*inv(R)*W*(H*A*x-am(iii));
    Lambda = Lambda  + A'*H'*W'*inv(R)*W*H*A;
    Omega = Omega  + A'*H'*W'*inv(R)*W*am(iii);
    if(tpm(jjj)<t(iii)-delay)
        dt = tpm(jjj)-t_state ;
        A = [1,dt,dt^2/2,0;0,1,dt,0;0,0,1,0;0,0,0,1];
        H = [1,0,0,0];
        R = rou_k_p^2;
        W =1;
        %E = E + (H*A*x-pm(iii))'*W'*inv(R)*W*(H*A*x-pm(iii));
        Lambda = Lambda  + A'*H'*W'*inv(R)*W*H*A;
        Omega = Omega  + A'*H'*W'*inv(R)*W*pm(jjj);
        jjj = jjj +1 ;
    end
    
    P = inv(Lambda);
    x = P*Omega;
    x_history(:,:,iii) = x;
    P_history(:,:,iii) = P;
    
end

figure,
title('Diff between gt to prediction at end point')
subplot(221)

plot(t,permute(x_history(1,1,:),[1,3,2])-p(end),'g');
grid on;
hold on;
std3 = permute(3*P_history(1,1,:).^0.5,[1,3,2]);
plot(t,std3,'r');
plot(t,-std3,'r');
legend('diff','3sdt');
ylim([-1,1]*median(std3)*3);
title('position diff');

subplot(222)

plot(t,permute(x_history(2,1,:),[1,3,2])-v(end),'g');
grid on;
hold on;
std3 = permute(3*P_history(2,2,:).^0.5,[1,3,2]);
plot(t,std3,'r');
plot(t,-std3,'r');
legend('diff','3sdt');
ylim([-1,1]*median(std3)*3);
title('velocity diff');


subplot(223)

plot(t,permute(x_history(3,1,:)+x_history(4,1,:),[1,3,2])-pdd-b_a,'g');
grid on;
hold on;
std3 = permute(3*P_history(3,3,:).^0.5,[1,3,2]);
plot(t,std3,'r');
plot(t,-std3,'r');
legend('diff','3sdt');
ylim([-1,1]*median(std3)*3);
title('acceleration diff');


subplot(224)

plot(t,permute(x_history(4,1,:),[1,3,2])-b_a,'g');
grid on;
hold on;
std3 = permute(3*P_history(4,4,:).^0.5,[1,3,2]);
plot(t,std3,'r');
plot(t,-std3,'r');
legend('diff','3sdt');
ylim([-1,1]*median(std3)*3);
title('bias diff');





P = inv(Lambda);
x = P*Omega;
gt = [p(end);v(end);pdd;b_a];
dx = x - gt;
mahal  = dx'*inv(P)*dx

disp('      GroundTruth |  state   |   diff');
disp(['pos   | ',num2str(gt(1)), '  |  ' ,num2str(x(1)), '  |  ' ,num2str(dx(1))]);
disp(['vel   | ',num2str(gt(2)), '  |  ' ,num2str(x(2)), '  |  ' ,num2str(dx(2))]);
disp(['acc   | ',num2str(gt(3)), '  |  ' ,num2str(x(3)), '  |  ' ,num2str(dx(3))]);
disp(['bias   | ',num2str(gt(4)), '  |  ' ,num2str(x(4)), '  |  ' ,num2str(dx(4))]);

disp(['mahalanobis distance = ' ,num2str(mahal)]);
