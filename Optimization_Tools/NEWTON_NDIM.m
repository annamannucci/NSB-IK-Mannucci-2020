%% NEWTON (N-DIM)
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
c1    = 0.0001;    % c1 > 0            
c2    = 0.01;      % c2 > 0
p     = 20;        % p >= 2
q     = 40;        % q >= 3
% ------------------------------------------------------------------------------------------------------ %

N     = size(X0,1);
% definition of the domain for N = 2 (necessary for the figures): "square domain"
% ------------------------------------------------------------------------------------------------------ %
xmin  = -6;             % minimum coordinate
xmax  =  6;             % maximum coordinate
dx    = 0.5;            % step between two consecutive coordinates
% ------------------------------------------------------------------------------------------------------ %

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
% If you want, you can add another inputs for the function armijo (line 71)

% gradient and hessian of objective function
F        = matlabFunction(Fsym,'vars',{x});
gradFsym = gradient(Fsym,x);
gradF    = matlabFunction(gradFsym,'vars',{x});
Qsym     = hessian(Fsym,x);  
Q        = matlabFunction(Qsym,'vars',{x});

X     = X0;               
Xit   = X';               
k     = 0;
                          
%% COMPUTATION
gradF_X0 = gradF(X0);
Q_X0     = Q(X0);

while (norm(gradF_X0)>tol)
               
          if (det(Q_X0)~=0)
              dN = -Q_X0\gradF_X0;
              condition1 = gradF_X0'*dN <= -c1*(norm(gradF_X0))^q;
              condition2 = (norm(dN))^p <= c2*norm(gradF_X0);
              if (condition1 && condition2)
                  d = dN;
              else
                  d = -gradF_X0;
              end
          else
              d = -gradF_X0;
          end

          % Chech the alpha value for the Armijo condition
          alfa_guess = 1;
          alfa = armijo(alfa_guess,X0,d,F,gradF);                    
      
          % UPDATE and STORAGE
          Xnew   = X0+alfa*d;
          X0     = Xnew;               
          Xitbis = [Xit;Xnew'];        
          Xit    = Xitbis;
          gradF_X0 = gradF(X0);
          Q_X0     = Q(X0);
          k     = k+1;
      
end

%% SOLUTION
Xmin = X0;              % point of minimum 
minf = F(X0);           % minimum value of the objective function

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