clear all;
close all;
clc;

%set the test you want to simulate: 
tests = ["SR"; "RP"; "iCAT-TCP"; "iCAT-SR"; "iCAT-conservative-RP"; "iCAT-RP"];

%"ST", Standard Approach with standard activation
%"iCAT-TCP", Standard Approach with regularization (cit. Simetti)
%"SR", Singularity Robust with standard activation matrix
%"iCAT-SR", Regularized Singularity Robust
%"RP", Reverse Priority with standard activation matrix
%"iCAT-conservative-RP",
%"iCAT-RP", Reverse Priority with regularization

%Directory to save data
directory = "/home/anna/Desktop/tests/obs10";

%Used just for naming
tasks = ["Collision Avoidance"; "Centroid"; "Perimeter"];

%enable plots and analysis
get_new_simulation = true;
plots = true;

%decide the set of obstacles to use
default_obs = false;
load('obstacles/obs10');
emax = zeros(length(tests),length(tasks)+1);
emean = zeros(length(tests),length(tasks)+1);
estd = zeros(length(tests),length(tasks)+1);

%% Start simulations
if get_new_simulation
 for t = 1 : length(tests)
 clearvars -except n_agents pobs default_obs dsafe_ths directory tests tasks epsilon_damp lambda_damp plots get_new_simulation k_vo t
 RTPflag = handle_name(tests(t));
 disp("Start simulating " + tests(t));

 entrapment
 save(directory+"/"+tests(t)+".mat", '-regexp', '^(?!(tests|tasks|t)$).');
 end
end

%% Plot Trajectories

figure('Name','Obstacle avoidance','Color',[1 1 1])
for t = 1 : length(tests)
 clearvars -except h_subplot n_agents n_circ dsafe_ths directory tests tasks epsilon_damp lambda_damp plots get_new_simulation t
 load(directory+"/"+tests(t)+".mat")
 if t == 7
   h_subplot(t) = subplot(3,3,t+2);
 else
  h_subplot(t) = subplot(3,3,t);
 end
 x = zeros(n_agents,length(q.time));
 y = zeros(n_agents,length(q.time));
 for i = 1 : 1 : n_agents
 for k = 1 : 1 : length(q.time)
 x(i,k) = q.signals.values(k,(i-1)*2+1);
 y(i,k) = q.signals.values(k,(i-1)*2+2);
 end
 plot(x(i,:), y(i,:),'.','LineWidth',1);
 hold on;
 end
 title(tests(t),'Interpreter','latex','FontSize',18)
 xlabel('x (m)','Interpreter','latex','FontSize',14);
 ylabel('y (m)','Interpreter','latex','FontSize',14);
 axis equal
 grid on
 for i = 1 : size(pobs,1)
 hold on
 viscircles([pobs(i,1),pobs(i,2)],dsafe);
 hold on;
 viscircles([pobs(i,1),pobs(i,2)],dsafe+dsafe_ths);
 end
end
linkaxes(h_subplot,'xy')
savefig(directory+"/trajectories.fig");
%% Merged velocities

figure('Name','Velocities','Color',[1 1 1])
for t = 1 : length(tests)
 clearvars -except h_subplot n_agents n_circ dsafe_ths directory tests tasks epsilon_damp lambda_damp plots get_new_simulation t
 load(directory+"/"+tests(t)+".mat")
 if t == 7
    h_subplot(t) = subplot(3,3,t+2);
 else
    h_subplot(t) = subplot(3,3,t);
 end
 plot(qdot.time, qdot.signals.values);
 title(tests(t),'Interpreter','latex','FontSize',16)
 xlabel('Time (sec)','Interpreter','latex','FontSize',14);
 ylabel('$\dot{q}$ ($m/s$)','Interpreter','latex','FontSize',14);
 ylim([-2 4])
 xlim([90 160])
 grid on
end
linkaxes(h_subplot,'xy')

savefig(directory+"/velocities.fig");
%% Results analysis

% plot obstacles
if (1)
     figure('Name', 'Collision Avoidance Static Obstacle');
     str = [];

     for t = 1 : length(tests)
     clearvars -except n_agents pobs default_obs dsafe_ths directory tests tasks epsilon_damp lambda_damp plots get_new_simulation t str
     load(directory+"/"+tests(t)+".mat")

     rel_pose = zeros(length(q.time),2);
     distance_ = zeros(length(q.time), n_agents*size(pobs,1));
     collision_ = zeros(length(q.time), n_agents*size(pobs,1));
     for k = 1 : length(q.time)
        for j = 1 : size(pobs,1)
            for i = 1 : n_agents
                rel_pose(k,:) = q.signals.values(k,(i-1)*2+1:(i-1)*2+2)-pobs(j,1:2); 
                 distance_(k,(j-1)*n_agents+i) = sqrt(rel_pose(k,1)^2+rel_pose(k,2)^2)-dsafe;
                 if (distance_(k,(j-1)*n_agents+i) > 0)
                    collision_(k,(j-1)*n_agents+i) = 0;
                 else
                    collision_(k,(j-1)*n_agents+i) = abs(distance_(k,(j-1)*n_agents+i));
                 end
            end
        end
     end
     collision_ = collision_';
     err = sum(collision_,1);
     str = [str;tests(t)];
     plot(q.time,err)
     hold on
     end
     plot([q.time(1) q.time(end)],[0 0],'r','linewidth',2)
     hold on
     str = [str;"desired value"];
     title('Collision avoidance task error ($m$)','Interpreter','latex');
     grid on
     xlabel('Time (sec)','Interpreter','latex');
     ylabel('$\sum\limits_{i=1,\dots,n}{\max{\left[0,\,d_{safe}-d_{i,obs_j}\right]}}$','Interpreter','latex');
     legend(str)
end

savefig(directory+"/task_collision_avoidance.fig");

clear rel_pose distance_ collision_ k j i err 
%% plot centroid error norm 

if plots
 figure('Name', 'Centroid Behaviour Error');
 str = [];

 for t = 1 : length(tests)
     clearvars -except n_agents pobs default_obs dsafe_ths directory tests tasks epsilon_damp lambda_damp plots get_new_simulation t str
     load(directory+"/"+tests(t)+".mat")
     err_center = qc_ref.signals.values(:,1:2)-center;
     err_norm = sqrt(sum(err_center.^2,2));
     emax(2) = max(err_norm(find(q.time == 90):end,:));
     emean(2) = mean(err_norm(find(q.time == 90):end,:));
     estd(2) = std(err_norm(find(q.time == 90):end,:));
     str = [str;tests(t)];
     plot(q.time,err_norm)
     hold on
     clear err_center err_norm
     save(directory+"/"+tests(t)+".mat", '-regexp', '^(?!(tests|tasks|t|str)$).');
 end
 grid on
 title('Task error centroid ($m$)','Interpreter','latex');
 ylabel('$\left|\left|q_{c,ref}-\frac{1}{n}\sum_{i=1}^{n}{q_i}\right|\right|$','Interpreter','latex');
 xlabel('Time (sec)','Interpreter','latex');
 legend(str)
end
savefig(directory+"/task_centroid.fig");

%% Perimeter
if plots

 figure('Name','Perimeter Behaviour','Color',[1 1 1])
 str = [];

 for t = 1 : length(tests)
 clearvars -except n_agents pobs default_obs dsafe_ths directory tests tasks epsilon_damp lambda_damp plots get_new_simulation t hAxB hAxS str
 load(directory+"/"+tests(t)+".mat");
 str = [str;tests(t)];

 perimeterq = zeros(length(q.time),1);
 error_perimeter = zeros(length(q.time),1);
 for k = 1 : length(q.time)
     perimeterq(k,:) = .5*(q.signals.values(k,1:2)'-q.signals.values(k,2*(n_agents-1)+1:2*(n_agents-1)+2)')'*(q.signals.values(k,1:2)'-q.signals.values(k,2*(n_agents-1)+1:2*(n_agents-1)+2)');
     for i = 1 : n_agents
         if (i < n_agents)
         perimeterq(k,:) = perimeterq(k,:) + .5*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-q.signals.values(k,2*i+1:2*i+2)')'*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-q.signals.values(k,2*i+1:2*i+2)');
         end
     end
     error_perimeter(k,:) = norm(perimeterq(k,:)-.5*n_agents*l^2); 
 end

 plot(q.time, error_perimeter)
 hold on
 emax(3) = max(error_perimeter(find(q.time == 90):end,:));
 emean(3) = mean(error_perimeter(find(q.time == 90):end,:));
 estd(3) = std(error_perimeter(find(q.time == 90):end,:));
 clear perimeterq error_perimeter k
 save(directory+"/"+tests(t)+".mat", '-regexp', '^(?!(tests|tasks|t|str)$).');
 end
 
grid on
title('Task error perimeter ($m^2$)','Interpreter','latex','FontSize',12)
ylabel('$\frac{1}{2}\left(d_{1,n}^2+\sum_{i=2}^n{d_{i,i-1}^2}-nl^2\right)$','Interpreter','latex','FontSize',12)
xlabel('Time (sec)','Interpreter','latex','FontSize',12);
legend(str,'Position',[0.913780664587004,0.798670511415543,0.082251081444742,0.123778498231006]);  
end

savefig(directory+"/task_perimeter.fig");

%% plot circular
 
if plots
 figure('Name', 'Circular Behaviour Error');
 str = [];

 for t = 1 : length(tests)
     clearvars -except n_agents pobs default_obs dsafe_ths directory tests tasks epsilon_damp lambda_damp plots get_new_simulation t ax1 ax2 ax3 str
     load(directory+"/"+tests(t)+".mat");
     err_circle = zeros(n_agents,length(q.time));
     r = l/(2*sin(pi/n_agents));
     for k = 1 : length(q.time)
         for i = 1 : n_agents
          err_circle(i,k) = r^2/2 - .5*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-center(k,:)')'*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-center(k,:)');
         end
     end
     err_circle_norm = sum(sqrt(err_circle.^2),1);

     emax(4) = max(err_circle_norm(:,find(q.time == 90):end));
     emean(4) = mean(err_circle_norm(:,find(q.time == 90):end));
     estd(4) = std(err_circle_norm(:,find(q.time == 90):end));

     str = [str;tests(t)];
     plot(q.time, err_circle_norm)
     hold on
     grid on
     title('Task error circular ($m^2$)','Interpreter','latex')
     ylabel('$\left|\left|\frac{1}{2}\left(r^2-d_{i,center}^2\right)\right|\right|$','Interpreter','latex','FontSize',12)
     xlabel('Time (sec)');
     legend(str);
     
     clear err_circle err_circle_norm k
     save(directory+"/"+tests(t)+".mat", '-regexp', '^(?!(tests|tasks|t|str)$).');
 end
 savefig(directory+"/costfun_circular.fig");
 
end