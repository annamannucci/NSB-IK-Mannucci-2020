%% Simulation parameters

Tc      = .05;           % sampling time
Tsim    = 200;          % simulation time
ff      = 1 - .5*Tc;    % forgetting factor

n_agents = 9;           % number of agents

% starts: initial positions (x,y) for each agent
% starts = [5 10, ...
%           5 0, ...
%           5 -10, ...
%           -5 15, ...
%           -5 10, ...
%           -5 5, ...
%           -5 -5, ...
%           -5 -10, ...
%           -5 -15];

%check another configuration for starts
theta = 2*pi/n_agents;
starts = zeros(n_agents*2,1);
for i = 1 : n_agents
     starts(2*i-1:2*i) = 10*[cos(theta*(i-1)) sin(theta*(i-1))];
end

% starts = [5 10, ...
%           4 0, ...
%           2 -10, ...
%           -1 15, ...
%           -5 3, ...
%           -1 5, ...
%           -5 -4, ...
%           -4 -9, ...
%           -5 -15];


% load predefined obstacles and collision avoidance parameters pobs(i,:) := pose
% (x,y) of the obstacle (generated just one time using the
% generate_obstacles function)
% if ~exist('obstacles','dir')
%     mkdir obstacles;
%     generate_obstacles;
% end
% 
% if ~exist('default_obs','var')
%     default_obs = true;
% end
% if (default_obs)
%     load('obstacles\obs0.mat');
% end

pobs = [
      140 2.5;
      145 -5;
      150 5;
      167.5 2.5;
      167.5 -2.5;
      170 5;
      170 -5
        ];

dsafe = 1;      % internal (safety) radius of the obstacles

% external radius of the obstacles: diameter of activation buffer 
% (distance from the border of the obstacle used for activating the task)
dsafe_ths = 1;  

%% set global parameters
epsilon_damp = 0.1;
lambda_damp = 0.1;
k_vo = 1;               % used in regularization to weight the contribution of tasks which are in transition

%% Task parameters: DISTANCE FROM CENTER and PERIMETER
% parameters for DISTANCE FROM CENTER task
% load circular and perimeter parameters
r = 10;                             % distance from center (r = 4 causes the collision, 10 is fine)
l = 2*r*cos(pi/2-pi/n_agents);      % chord for minimizing the escape distances
tol = n_agents*r^2/50;

%% Task parameters: CENTROID
% Target reference trajectory for centroid motion
xa = [0 0];
xb = [200 0];
Tab = 180;
time = 0 : Tc : Tsim;
x = zeros(length(time),2);
xd = zeros(length(time),2);
xdd = zeros(length(time),2);

% We use a fifth order interpolating polynomial law (check Section IV
% in the paper: "On Null Space-Based Inverse Kinematics Techniques for 
% Fleet Management: Toward Time-Varying Task Activation" (Mannucci, A.,
% Caporale, D., Pallottino, L.).

M = [Tab^5 Tab^4 Tab^3;
    5*Tab^4 4*Tab^3 3*Tab^2;
    20*Tab^3 12*Tab^2 6*Tab];
b = [xb(1); 0; 0];
a = M\b;

clear M b

for i = 1 : length(time)
    if time(i) <= Tab
        x(i,1) = a(1)*time(i)^5+a(2)*time(i)^4+a(3)*time(i)^3;
        xd(i,1) = 5*a(1)*time(i)^4+4*a(2)*time(i)^3+3*a(3)*time(i)^2;
        xdd(i,1) = 20*a(1)*time(i)^3+12*a(2)*time(i)^2+6*a(3)*time(i);
    else
        x(i,1) = x(Tab/Tc+1,1,1);
        xd(i,1) = xd(Tab/Tc+1,1);
        xdd(i,1) = xdd(Tab/Tc+1,1);
    end
end

qc_ref.time = time';
qc_ref.dimensions = 2*n_agents;

qcdot_ref.time = time';
qcdot_ref.dimensions = 2*n_agents;

for index = 1 : length(qc_ref.time)
    qc_ref.signals.values(index,:) = x(index,:);
    qcdot_ref.signals.values(index,:) = xd(index,:);
end

clear xa xb x xd xdd Tab time index


%% Start simulation
tic
sim('entrapment_multi');
toc