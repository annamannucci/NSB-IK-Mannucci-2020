%% RANDOM WALK (N-DIM)
% ATTENTION: you can set the variables that appear into the interval
% between % --------- % !!!

clear all
close all
clc
format long

%% DATA
% ------------------------------------------------------------------------------------------------------ %
X0         = [5;5];     % coordinates of the first attempt (column vector)
tol        = 1e-4;      % tolerance to stop the algorithm                   (***advisable 1e-4/1e-5***)
% ------------------------------------------------------------------------------------------------------ %

N     = size(X0,1);
% definition of the domain for N = 2 (necessary for the figures): "square domain"
% ------------------------------------------------------------------------------------------------------ %
xmin  = -6;             % minimum coordinate
xmax  =  6;             % maximum coordinate
dx    = 0.5;            % step between two consecutive coordinates
% ------------------------------------------------------------------------------------------------------ %

%% INITIALIZATION
alfa0 = 1;                
% ------------------------------------------------------------------------------------------------------ %
ndir  = 2*10^(N-1);     % max number of directions assessed for each point  (***advisable 2*10^(N-1)***)
% ------------------------------------------------------------------------------------------------------ %
X     = X0;               
Xit   = X';               
alfa  = alfa0;            
k     = 0;                
                          
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

F = matlabFunction(Fsym,'vars',{x});

%% COMPUTATION
while (alfa>tol)
               
      r    = random('unif',-1,1,N,1);   
      d    = r./norm(r);                
      Xnew = X0+alfa*d;                 
      fnew = F(Xnew);    
      f0   = F(X0);      
      
      if fnew<f0                  
         X0     = Xnew;               
         k      = 0;                 
         Xitbis = [Xit;Xnew'];        
         Xit    = Xitbis;
      else
         k = k+1;      
      end
      
      if k>ndir               
         alfa = alfa*0.5;    
         k = 0;             
      end

end

%% SOLUTION
Xmin = X0;              % point of minimum 
minf = F(X0);           % minimum value of the objective function

%% PLOT

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