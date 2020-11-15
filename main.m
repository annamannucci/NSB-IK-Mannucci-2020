clear all;
close all;
clc;

%ATTENTION: this simulation is customized with this fixed setup of just 4
%tasks enabled. The simulation can be customized for different set of tasks
%by changing the NSB_IK_controller as follows:
%1) line 4 (number of tasks) -- n_tasks
%2) line 13 (number of subtasks) -- h_length_
%3) lines 92-104 to use the gradient for defining the circular task (see
%   Section IV.C)
%4) lines 143-176 (piling up the correct Jacobians)

%Set the test you want to simulate: 
%"ST", Standard Approach with standard activation
%"iCAT-TCP", Standard Approach with regularization (cit. Simetti)
%"SR", Singularity Robust with standard activation matrix
%"iCAT-SR", Regularized Singularity Robust
%"RP", Reverse Priority with standard activation matrix
%"iCAT-conservative-RP",
%"iCAT-RP", Reverse Priority with regularization
tests = ["ST"; "SR"; "RP"; "iCAT-TPC"; "iCAT-SR"; "iCAT-conservative-RP"; "iCAT-RP"];

%Directory to save data
directory = "tests/";

%enable simulation, plots and analysis
get_new_simulation = true;
plots = true;

%NOTE: parameters related to the simulation (e.g., obstacles, damping parameters,
%number of robots ... ) are set inside the entrapment file.
%% Start simulations
if get_new_simulation
 for t = 1 : length(tests)
 clearvars -except n_agents dsafe_ths directory tests epsilon_damp lambda_damp plots get_new_simulation k_vo t
 RTPflag = handle_name(tests(t));
 disp("Start simulating " + tests(t));
 entrapment
 save(directory+"/"+tests(t)+".mat", '-regexp', '^(?!(t|i)$).');
 end
end

%% Plot paths

for t = 1 : length(tests)
 clearvars -except h_subplot n_agents n_circ dsafe_ths directory tests epsilon_damp lambda_damp plots get_new_simulation t
 load(directory+"/"+tests(t)+".mat")
 figure('Name','Trajectories','Color',[1 1 1])

 %Plot obstacles
 for i = 1 : size(pobs,1)
   hold on
    viscircles([pobs(i,1),pobs(i,2)],dsafe);
    hold on;
    viscircles([pobs(i,1),pobs(i,2)],dsafe+dsafe_ths);
 end

 %Plot paths
 x = zeros(n_agents,length(q.time));
 y = zeros(n_agents,length(q.time));
 for i = 1 : 1 : n_agents
 for k = 1 : 1 : length(q.time)
    x(i,k) = q.signals.values(k,2*i-1);
    y(i,k) = q.signals.values(k,2*i);
 end
 plot(x(i,:), y(i,:),'.','LineWidth',1);
 hold on;
 end
 title(tests(t),'Interpreter','latex','FontSize',18)
 xlabel('x (m)','Interpreter','latex','FontSize',14);
 ylabel('y (m)','Interpreter','latex','FontSize',14);
 axis equal
 grid on
end

%% Plot velocities
for t = 1 : length(tests)
 clearvars -except h_subplot n_agents n_circ dsafe_ths directory tests epsilon_damp lambda_damp plots get_new_simulation t
 load(directory+"/"+tests(t)+".mat")
 figure('Name','Velocities','Color',[1 1 1])
 plot(qdot.time, qdot.signals.values);
 title(tests(t),'Interpreter','latex','FontSize',16)
 xlabel('Time (sec)','Interpreter','latex','FontSize',14);
 ylabel('$\dot{q}$ ($m/s$)','Interpreter','latex','FontSize',14);
 ylim([-2 4])
 grid on
end

%% Results analysis

% plot obstacles
if plots
 figure('Name', 'Collision Avoidance Static Obstacle');
 str = [];

 for t = 1 : length(tests)
 clearvars -except n_agents dsafe_ths directory tests epsilon_damp lambda_damp plots get_new_simulation t str
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
%% plot centroid error norm 

if plots
 figure('Name', 'Centroid Behaviour Error');
 str = [];

 for t = 1 : length(tests)
     clearvars -except n_agents dsafe_ths directory tests epsilon_damp lambda_damp plots get_new_simulation t str
     load(directory+"/"+tests(t)+".mat")
     err_center = qc_ref.signals.values(:,1:2)-center;
     err_norm = sqrt(sum(err_center.^2,2));
     max(err_norm(find(q.time == 90):end,:))
     mean(err_norm(find(q.time == 90):end,:))
     std(err_norm(find(q.time == 90):end,:))
     str = [str;tests(t)];
     plot(q.time,err_norm)
     hold on
     save(directory+"/"+tests(t)+".mat",'err_norm','-append');
 end


 grid on
 title('Task error centroid ($m$)','Interpreter','latex');
 ylabel('$\left|\left|q_{c,ref}-\frac{1}{n}\sum_{i=1}^{n}{q_i}\right|\right|$','Interpreter','latex');
 xlabel('Time (sec)','Interpreter','latex');
 legend(str)
end

%% plot circular
 
if plots
 figure('Name', 'Circular Behaviour Error');
 ax1 = subplot(2,2,[1,2]);
 ax2 = subplot(2,2,3);
 ax3 = subplot(2,2,4);
 str = [];

 for t = 1 : length(tests)
     clearvars -except n_agents dsafe_ths directory tests epsilon_damp lambda_damp plots get_new_simulation t ax1 ax2 ax3 str
     load(directory+"/"+tests(t)+".mat");
     err_circle = zeros(n_agents,length(q.time));
     err_circle_distance = zeros(n_agents,length(q.time));
     err_circle_mean_distance = zeros(length(q.time),1);
     err_circle_max_distance = zeros(length(q.time),1);
     r = l/(2*sin(pi/n_agents));
     for k = 1 : length(q.time)
         for i = 1 : n_agents
          err_circle(i,k) = r^2/2 - .5*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-center(k,:)')'*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-center(k,:)');
          err_circle_distance(i,k) = abs(r - norm(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-center(k,:)'));
         end
         err_circle_mean_distance(k) = 1/n_agents*sum(err_circle_distance(:,k));
         err_circle_max_distance(k) = max(err_circle_distance(:,k));
     end
     err_circle_norm = sum(sqrt(err_circle.^2),1);

     %max(err_circle_norm(:,find(q.time == 90):end))
     %mean(err_circle_norm(:,find(q.time == 90):end))
     %std(err_circle_norm(:,find(q.time == 90):end))

     str = [str;tests(t)];
     hold(ax1,'on');
     hold(ax2,'on');
     hold(ax3,'on');
     plot(ax1, q.time, err_circle_norm)
     plot(ax2, q.time, err_circle_mean_distance)
     plot(ax3, q.time, err_circle_max_distance)

     grid(ax1,'on');
     title(ax1,'Task error circular ($m^2$)','Interpreter','latex')
     ylabel(ax1,'$\left|\left|\frac{1}{2}\left(r^2-d_{i,center}^2\right)\right|\right|$','Interpreter','latex','FontSize',12)
     xlabel(ax1,'Time (sec)');
     legend(ax1,str);

     grid(ax2,'on');
     title(ax2,'Mean radial error between agents (m)','Interpreter','latex');
     ylabel(ax2,'$\frac{1}{n}\sum_{i=1}^{n}{\left|d_{i,center}-r\right|}$','Interpreter','latex','FontSize',12)
     xlabel(ax2,'Time (sec)');
     legend(ax2,str);

     grid(ax3,'on');
     title(ax3,'Max radial error between agents (m)','Interpreter','latex');
     ylabel(ax3,'$\max_{i=1,\dots ,n}{\left|d_{i,center}-r\right|}$','Interpreter','latex','FontSize',12)
     xlabel(ax3,'Time (sec)');
     legend(ax3,str);
     save(directory+"/"+tests(t)+".mat",'err_circle_norm','-append');
 end
 

 clearvars ax1 ax2 ax3
 
end

%% Perimeter
if plots

 figure('Name','Perimeter Behaviour','Color',[1 1 1])
 %# build axes positions
 hBig = [subplot(2,2,[1 2]) subplot(2,2,3) subplot(2,2,4)]; %# create subplots
 posBig = get(hBig, 'Position'); %# record their positions
 delete(hBig) %# delete them
 posSmall{1} = [0.45, 0.65, 0.4, 0.25];

 %# create axes (big/small)
 hAxB(1) = axes('Position',posBig{1});
 hAxB(2) = axes('Position',posBig{2});
 hAxB(3) = axes('Position',posBig{3});
 hAxS(1) = axes('Position',posSmall{1});
 str = [];

 for t = 1 : length(tests)
 clearvars -except n_agents dsafe_ths directory tests epsilon_damp lambda_damp plots get_new_simulation t hAxB hAxS str
 load(directory+"/"+tests(t)+".mat");
 str = [str;tests(t)];

 perimeterq = zeros(length(q.time),1);
 error_perimeter = zeros(length(q.time),1);
 error_sides = zeros(length(q.time),1);
 error_max = zeros(length(q.time),1);
 for k = 1 : length(q.time)
     perimeterq(k,:) = .5*(q.signals.values(k,1:2)'-q.signals.values(k,2*(n_agents-1)+1:2*(n_agents-1)+2)')'*(q.signals.values(k,1:2)'-q.signals.values(k,2*(n_agents-1)+1:2*(n_agents-1)+2)');
     error_sides(k,:) = abs(sqrt((q.signals.values(k,1:2)'-q.signals.values(k,2*(n_agents-1)+1:2*(n_agents-1)+2)')'*(q.signals.values(k,1:2)'-q.signals.values(k,2*(n_agents-1)+1:2*(n_agents-1)+2)'))-l);
     error_max(k,:) = abs(sqrt((q.signals.values(k,1:2)'-q.signals.values(k,2*(n_agents-1)+1:2*(n_agents-1)+2)')'*(q.signals.values(k,1:2)'-q.signals.values(k,2*(n_agents-1)+1:2*(n_agents-1)+2)'))-l);
     for i = 1 : n_agents
         if (i < n_agents)
         perimeterq(k,:) = perimeterq(k,:) + .5*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-q.signals.values(k,2*i+1:2*i+2)')'*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-q.signals.values(k,2*i+1:2*i+2)');
         error_sides(k,:) = error_sides(k,:) + abs(sqrt((q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-q.signals.values(k,2*i+1:2*i+2)')'*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-q.signals.values(k,2*i+1:2*i+2)'))-l);
         tmp = abs(sqrt((q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-q.signals.values(k,2*i+1:2*i+2)')'*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-q.signals.values(k,2*i+1:2*i+2)'))-l);
             if (tmp > error_max(k,:))
             error_max(k,:) = tmp;
             end
         end
     end
     error_sides(k,:) = error_sides(k,:)/n_agents;
     error_perimeter(k,:) = norm(perimeterq(k,:)-.5*n_agents*l^2); 
 end

 hold(hAxB(1),'on')
 hold(hAxB(2),'on')
 hold(hAxB(3),'on')
 hold(hAxS(1),'on')
 plot(hAxB(1),q.time, error_perimeter)
 indexOfInterest = (q.time < 200) & (q.time > 100); % range of t near perturbation
 plot(hAxS(1),q.time(indexOfInterest),error_perimeter(indexOfInterest)) % plot on new axes

 plot(hAxB(2),q.time, error_sides)
 plot(hAxB(3),q.time, error_max)
 
 %max(error_perimeter(find(q.time == 90):end,:))
 %mean(error_perimeter(find(q.time == 90):end,:))
 %std(error_perimeter(find(q.time == 90):end,:))
 save(directory+"/"+tests(t)+".mat",'error_perimeter','-append');
 end
 
 grid(hAxS(1),'on');
 
 grid(hAxB(1),'on');
 title(hAxB(1),'Task error perimeter ($m^2$)','Interpreter','latex','FontSize',12)
 ylabel(hAxB(1),'$\frac{1}{2}\left(d_{1,n}^2+\sum_{i=2}^n{d_{i,i-1}^2}-nl^2\right)$','Interpreter','latex','FontSize',12)
 xlabel(hAxB(1),'Time (sec)','Interpreter','latex','FontSize',12);
 legend(hAxB(1),str,'Position',[0.913780664587004,0.798670511415543,0.082251081444742,0.123778498231006]); 
 
 grid(hAxB(2),'on');
 title(hAxB(2),'Mean error for each side between agents (m)','Interpreter','latex','FontSize',12);
 ylabel(hAxB(2),'$\frac{1}{n}\left(\left|d_{1,n}-l\right|+\sum_{i=2}^n{\left|d_{i,i-1}-l\right|}\right)$','Interpreter','latex','FontSize',12)
 xlabel(hAxB(2),'Time (sec)','Interpreter','latex','FontSize',12);
 
 grid(hAxB(3),'on');
 title(hAxB(3),'Max error for each side between agents (m)','Interpreter','latex','FontSize',12);
 ylabel(hAxB(3),'$\max_{i=1,\dots ,n}{\left|d_{i,i-1}-l\right|}$','Interpreter','latex','FontSize',12)
 xlabel(hAxB(3),'Time (sec)','Interpreter','latex','FontSize',12);
 
 clearvars hAxB hAxS
 
end

 

