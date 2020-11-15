%% plot centroid error norm 

clearvars -except n_agents dsafe_ths directory teRPs tasks epsilon_damp lambda_damp plots get_new_simulation t RPr
load(directory+"/RP.mat")
err_center_1 = qc_ref.signals.values(:,1:2)-center;
err_norm_1 = sqrt(sum(err_center_1.^2,2));
min(err_norm_1)

clearvars -except err_norm_1 n_agents dsafe_ths directory teRPs tasks epsilon_damp lambda_damp plots get_new_simulation t RPr
load(directory+"/RP1.mat")
err_center_2 = qc_ref.signals.values(:,1:2)-center;
err_norm_2 = sqrt(sum(err_center_2.^2,2));
min(err_norm_2)

ratio = zeros(length(err_norm_2),1);
for i = 1 : length(err_norm_2)
    ratio(i) = min(err_norm_1(i)/1e-3,10);
    if ratio(i) < 1
        break
    end
end
q.time(i)


%% plot circular
 clearvars -except n_agents dsafe_ths directory teRPs tasks epsilon_damp lambda_damp plots get_new_simulation t ax1 ax2 ax3 RPr
 load(directory+"/RP.mat");
 err_circle = zeros(n_agents,length(q.time));
 r = l/(2*sin(pi/n_agents));
 for k = 1 : length(q.time)
     for i = 1 : n_agents
      err_circle(i,k) = r^2/2 - .5*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-center(k,:)')'*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-center(k,:)');
     end
 end
 err_circle_norm_1= sum(sqrt(err_circle.^2),1);
 min(err_circle_norm_1)

 load(directory+"/RP1.mat");
 err_circle = zeros(n_agents,length(q.time));
 r = l/(2*sin(pi/n_agents));
 for k = 1 : length(q.time)
     for i = 1 : n_agents
      err_circle(i,k) = r^2/2 - .5*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-center(k,:)')'*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-center(k,:)');
     end
 end
 err_circle_norm_2= sum(sqrt(err_circle.^2),1);
 min(err_circle_norm_2)
 
ratio = zeros(length(err_circle_norm_2),1);
for i = 1 : length(err_circle_norm_2)
    ratio(i) = min(err_circle_norm_1(i)/err_circle_norm_2(i),10);
    if ratio(i) < 2
        break
    end
end
q.time(i)

%% Perimeter
 clearvars -except n_agents dsafe_ths directory teRPs tasks epsilon_damp lambda_damp plots get_new_simulation t hAxB hAxS RPr
 load(directory+"/RP.mat");

 perimeterq = zeros(length(q.time),1);
 error_perimeter_1 = zeros(length(q.time),1);
 for k = 1 : length(q.time)
     perimeterq(k,:) = .5*(q.signals.values(k,1:2)'-q.signals.values(k,2*(n_agents-1)+1:2*(n_agents-1)+2)')'*(q.signals.values(k,1:2)'-q.signals.values(k,2*(n_agents-1)+1:2*(n_agents-1)+2)');
     for i = 1 : n_agents
         if (i < n_agents)
         perimeterq(k,:) = perimeterq(k,:) + .5*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-q.signals.values(k,2*i+1:2*i+2)')'*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-q.signals.values(k,2*i+1:2*i+2)');
         end
     end
     error_perimeter_1(k,:) = norm(perimeterq(k,:)-.5*n_agents*l^2); 
 end
 min(error_perimeter_1)
 
 load(directory+"/RP1.mat");
 perimeterq = zeros(length(q.time),1);
 error_perimeter_2 = zeros(length(q.time),1);
 for k = 1 : length(q.time)
     perimeterq(k,:) = .5*(q.signals.values(k,1:2)'-q.signals.values(k,2*(n_agents-1)+1:2*(n_agents-1)+2)')'*(q.signals.values(k,1:2)'-q.signals.values(k,2*(n_agents-1)+1:2*(n_agents-1)+2)');
     for i = 1 : n_agents
         if (i < n_agents)
         perimeterq(k,:) = perimeterq(k,:) + .5*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-q.signals.values(k,2*i+1:2*i+2)')'*(q.signals.values(k,2*(i-1)+1:2*(i-1)+2)'-q.signals.values(k,2*i+1:2*i+2)');
         end
     end
     error_perimeter_2(k,:) = norm(perimeterq(k,:)-.5*n_agents*l^2); 
 end
min(error_perimeter_2)

ratio = zeros(length(error_perimeter_1),1);
for i = 1 : length(error_perimeter_1)
    %ratio(i) = min(error_perimeter_1(i)/(sqrt(.5*n_agents)*1e-3),10);
    if error_perimeter_1(i) < sqrt(.5*n_agents)*1e-3
        break
    end
end
q.time(i)
