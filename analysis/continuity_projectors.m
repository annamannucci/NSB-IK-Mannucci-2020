close all
clear all
clc

ta = pi/6;
tb = pi/4;
X = [cos(ta) sin(ta); cos(tb) sin(tb)];
prA_active_B_disactive = pinv(X(1,:));
final = pinv(X);
prA_active_B_active = final(:,1);
prB_active_A_active = final(:,2);

%Fissato tb, come varia il proiettore al variare del fattore di attivazione, 
%considerando la dumped pinv o la pinvsc
prB_active_A_disactive = pinv(X(2,:));
eps = logspace(-12,0,10000);
N  = length(eps);
prA_Binv = zeros(2,N);
prB_Binv = zeros(2,N);
d = zeros(1,N);

epsilon = 1;
lambda = 1;

for i = 1 : N
    A = diag([1 eps(i)]);
    [M, d(i)] = damped_pinv(X,epsilon,lambda);
    prA_Binv(:,i) = M(:,1);
    prB_Binv(:,i) = M(:,2);
end

figure('Name','Projector 1')
subplot(2,1,1)
plot(0,prA_active_B_disactive(1),'o','Color','b')
hold on
plot(1,prA_active_B_active(1),'*','Color','b')
hold on
plot(eps,prA_Binv(1,:),'b')
hold on
plot(0,prA_active_B_disactive(2),'o','Color','r')
hold on
plot(1,prA_active_B_active(2),'*','Color','r')
hold on
plot(eps,prA_Binv(2,:),'r')
hold on
plot(0,atan2(prA_active_B_disactive(2),prA_active_B_disactive(1)),'*','Color','k')
hold on
plot(1,atan2(prA_active_B_active(2),prA_active_B_active(1)),'*','Color','k')
hold on
plot(eps,atan2(prA_Binv(2,:),prA_Binv(1,:)),'k')
%legend('x','y','atan(y,x)')

subplot(2,1,2)
plot(0,prB_active_A_disactive(1),'o','Color','b')
hold on
plot(1,prB_active_A_active(1),'*','Color','b')
hold on
plot(eps,prB_Binv(1,:),'b')
hold on
plot(0,prB_active_A_disactive(2),'o','Color','r')
hold on
plot(1,prB_active_A_active(2),'*','Color','r')
hold on
plot(eps,prB_Binv(2,:),'r')
hold on
plot(0,atan2(prB_active_A_disactive(2),prB_active_A_disactive(1)),'*','Color','k')
hold on
plot(1,atan2(prB_active_A_active(2),prB_active_A_active(1)),'*','Color','k')
hold on
plot(eps,atan2(prB_Binv(2,:),prB_Binv(1,:)),'k')