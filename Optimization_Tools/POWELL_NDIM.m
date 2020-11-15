%% POWELL (N-DIM)
% ATTENTION: you can set the variables that appear into the interval
% between % --------- % !!!

clear all
close all
clc
format long

%% DATA
% ------------------------------------------------------------------------------------------------------ %
X0    = [5;5];     % coordinates of the first attempt (column vector)
tol   = 1e-4;      % tolerance to stop the algorithm                   (***advisable 1e-4/1e-5***)
% ------------------------------------------------------------------------------------------------------ %

% definition of the domain for N = 2 (necessary for the figures): "square domain"
% ------------------------------------------------------------------------------------------------------ %
xmin  = -6;             % minimum coordinate
xmax  =  6;             % maximum coordinate
dx    = 0.5;            % step between two consecutive coordinates
% ------------------------------------------------------------------------------------------------------ %

% initialization
N      = size(X0,1);      % state's dimension
k      = 0;               % iterations' counter
check1 = 1;               % block condition for points' values
check2 = 1;               % block condition for function's values
Xk     = X0;              % current point
Xit    = X0';             % storage the points
alpha_guess = 1;          % guess value for Armijo-backtracking
Y      = zeros(N,N+1);    % set of points assessed for each cycle
reduct = zeros(N,1);      % vector of the progressive reductions of the objective function
S      = eye(N);          % main directions (column vector) at the first cycle

%% OBJECTIVE FUNCTION_HANDLE
x    = sym('x',[N 1]);
% NOTE: before to continue, you have to insert the objective function,
% depending by variables x(1), x(2), etc.
% Example 1: Fsym = (x(1)-1).^4+x(2).^2.*(x(1)-2).^2+(x(3)+1).^2;
% Example 2: Fsym = (x(1)-2).^4+x(2).^2.*(x(1)-2).^2+(x(2)+1).^2;     
% ------------------------------------------------------------------------------------------------------ %
Fsym = (x(1)-2).^4+x(2).^2.*(x(1)-2).^2+(x(2)+1).^2;     
% ------------------------------------------------------------------------------------------------------ %

%% STOP FOR INPUT DATA
% If you want, you can add another inputs for the function backtr (line 77)

F        = matlabFunction(Fsym,'vars',{x});

%% COMPUTATION

while (check1 || check2)
    
    % DIRECTION with MAX REDUCTION
    Y(:,1) = Xk;
    for j = 1:N
        sj  = S(:,j);   % j   <-- direction
        yj_ = Y(:,j);   % j-1 <-- point
        alphaj = backtr(alpha_guess,yj_,sj,F);
        Y(:,j+1)  = Y(:,j)+alphaj.*sj;
        if j<N
            reduct(j) = F(Y(:,j))-F(Y(:,j+1));  
        end
    end
    % max reduction
    [DELTA,index_DELTA] = max(reduct);

    F1 = F(Y(:,1));
    F2 = F(Y(:,end));
    F3 = F(2.*Y(:,end)-Y(:,1));
    
    % NEW POINT
    if ((F3>=F1) || (((F1-2*F2+F3)*(F1-F2-DELTA)^2)>=0.5*DELTA*(F1-F3)^2))
        % sj(k+1)==sj(k) --> same directions
        % X(k+1)==Y(N)   --> to start from Y(:,end)
        Xk_new = Y(:,end);
    else
        dk = (Y(:,end)-Xk)/norm(Y(:,end)-Xk);
        alphaj = backtr(alpha_guess,Xk,dk,F);
        Xk_new = Xk+alphaj.*dk;
        % delete the direction index_DELTA
        S(:,index_DELTA) = [];
        % input the new direction dk
        Snew = [S,dk];
        S = Snew;  
    end
    
    % UPDATE and STORAGE
    k      = k+1;
    Xitbis = [Xit;Xk_new'];
    Xit    = Xitbis;
    
    % CHECK
    check1 = (norm(Xk_new-Xk)/max(1,norm(Xk)))>tol;
    check2 = (abs(F(Xk_new)-F(Xk))/(max(1,norm(Xk))))>tol;
    
    Xk = Xk_new;
end

Xmin  = Xk;
iter  = k;
Fmin  = F(Xk);
punti = size(Xit,1);

%% PLOT
clc

if N==2
    
    dom   = xmin:dx:xmax;        
    Fbis  = matlabFunction(Fsym,'vars',{x(1),x(2)});
    [x,y] = meshgrid(dom);       
    z = Fbis(x,y);

    % 1
    figure('Name','Function 3D + Level Lines');
    meshc(x,y,z);
    xlabel(' x ');
    ylabel(' y ');
    zlabel(' f(x,y) ');
    grid on;
    axis square;
    title('Function 3D + Level Lines');

    % 2 - option 1
    punti = size(Xit,1);
    figure('Name','Level Lines');
    [C,H] = contour(x,y,z,10);
    hold on;
    colorbar;
    plot(Xit(:,1),Xit(:,2),'o--k','LineWidth',1.5,'MarkerSize',6);
    grid on;
    axis on;
    title(sprintf('Points assessed:%d - [Xmin=%g] - [Ymin=%g]',punti,Xmin(1),Xmin(2)));
    clabel(C,H);

    % 2 - option 2
    [gradx,grady] = gradient(z,dx);
    figure('Name','Level Lines');
    contourf(x,y,z,20);
    hold on;
    colorbar;
    plot(Xit(:,1),Xit(:,2),'o--k','LineWidth',1.5,'MarkerSize',6,'MarkerFaceColor','k');
    quiver(x,y,gradx,grady,'b','LineWidth',1);
    axis on;
    title(sprintf('Points assessed:%d - [Xmin=%g] - [Ymin=%g]',punti,Xmin(1),Xmin(2)));
    xlabel('x');
    ylabel('y');

end

figure('Name','Points')
for j=1:N
    subplot(N,1,j);
    plot(0:size(Xit,1)-1,Xit(:,j),'o-b','LineWidth',1.5,'MarkerSize',6);
    set(gca,'XMinorTick','on','YMinorTick','on');
    ylabel(sprintf('component-x%d',j));
    legend(sprintf('x%d_m_i_n=%g',j,Xmin(j)));
    grid on;
    axis on;
    if j==N
        xlabel(' iterations [int] ');
    end
end