n_obs = size(pobs,1);
J1 = zeros(n_obs*n_agents,2*n_agents); %preallocation of Jacobians
xdot1 = zeros(n_obs*n_agents,1);       %preallocation of tasks

for j = 1 : n_obs
    for i = 1 : 1 : n_agents
        %Task: distance from the obstacle greater than a threshold
        dcurr = sqrt((q(2*i-1,1)-pobs(j,1))^2 + (q(2*i,1)-pobs(j,2))^2);
        J1((j-1)*n_agents+i,2*i-1:2*i) = (1/dcurr)*[q(2*i-1,1)-pobs(j,1) q(2*i,1)-pobs(j,2)]; %distance between obstacle and robot the goal
        xdot1((j-1)*n_agents+i,1) = .8 * (dsafe+dsafe_ths-dcurr);
    end
end


J6 = zeros(2,2*n_agents);
for i = 1 : n_agents
    J6(:,(i-1)*2+1:2*i) = 1/n_agents*eye(2);                   
end
%se si rimuove la velocità non è garantita la convergenza a zero dell'errore
xdot6 = n_agents*xdot6;
J6 = n_agents*J6;


%% 2 agenti
n_agents = 3;
J6 = zeros(2,2*n_agents);
for i = 1 : n_agents
    J6(:,(i-1)*2+1:2*i) = 1/n_agents*eye(2);                   
end

th = sym('th',[2*n_agents,1],'real');
J1 = [cos(th(1)) sin(th(1)) 0 0;
      cos(th(2)) sin(th(2)) 0 0];
syms h1 h2 real
H1 = diag([h1,h2]);
simplify(svd(J1))
svd(J6'*J6) %1/2 e 1/2;
svd(J1'*H1*J1+J6'*J6)
simplify(svd(J1'*H1*J1))

                                                                                                                                                                             0
 (2^(1/2)*(h1*h2 - h1*(h1^2 + 2*cos(2*th1 - 2*th2)*h1*h2 + h2^2)^(1/2) - h2*(h1^2 + 2*cos(2*th1 - 2*th2)*h1*h2 + h2^2)^(1/2) + h1^2 + h2^2 + h1*h2*cos(2*th1 - 2*th2))^(1/2))/2
 (2^(1/2)*(h1*h2 + h1*(h1^2 + 2*cos(2*th1 - 2*th2)*h1*h2 + h2^2)^(1/2) + h2*(h1^2 + 2*cos(2*th1 - 2*th2)*h1*h2 + h2^2)^(1/2) + h1^2 + h2^2 + h1*h2*cos(2*th1 - 2*th2))^(1/2))/2
 
J1 = [cos(pi/4) sin(pi/4) 0 0 0 0;
      cos(pi/8) sin(pi/8) 0 0 0 0;
      0 0 cos(pi/3) sin(pi/3) 0 0;
      0 0 cos(pi) sin(pi) 0 0];
  svd(J1)
  svd([J1;J6])
  
  J1 = [cos(th(1)) sin(th(1));
      cos(th(2)) sin(th(2))];
  
  alpha = linspace(-pi,pi,360);
  sp = zeros(length(alpha),2);
  for i = 1:length(alpha)
      sp(i,1) = (2^(1/2)*(2 - 2^(1/2)*(cos(2*alpha(i)) + 1)^(1/2))^(1/2))/2;
      sp(i,2) = (2^(1/2)*(2^(1/2)*(cos(2*alpha(i)) + 1)^(1/2) + 2)^(1/2))/2;
  end
  plot(th,sp(:,1),th,sp(:,2));
%%