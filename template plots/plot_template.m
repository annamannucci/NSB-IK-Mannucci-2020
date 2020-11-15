%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% TEMPLATE FOR NICE PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suggestions
% - reduce the sampling of your data to reduce the final file size

%% Load data
load('./../data/data_test')

%% General variables

save_ck=1; % flag for saving the plots
% colors
b = [0.0000    0.4470    0.7410];
r = [0.8500    0.3250    0.0980];
y = [0.9290    0.6940    0.1250];
p = [0.4940    0.1840    0.5560];
p = [0.25 0.25 0.25];
g = [0.4660    0.6740    0.1880];
lb = [0.3010    0.7450    0.9330];  % light blue
dr = [0.6350    0.0780    0.1840];  % dark red

% Plot parameters
LineWidth = 1;
FontSize = 18;
FontSize_label = 14;
FontSize_axes = 12;
MarkerSize = 6;
Position = [800 60];        % position of the figure on the screen
Dimension = [880 600];      % deimension of the figure on the screen
xTick_step = 1; 

dim_raw = 3;
dim_col = 2;
dim_fig = [dim_raw, dim_col];
marg_axes = [.01 .079];     % margin between plots: Heigth, Width
marg_hBox = [.07 .01];      % margin of the big box: Top, Botton
marg_wBox = [.085 .017];    % margin of the big box: Left, Right  

% useful variables 
r2d = 180/pi;
d2r = pi/180;

% variables to allign the labels along the y axis
ylabh = zeros(dim_raw, dim_col);
ylabelPosX = zeros(dim_raw, dim_col,3);

% Y margin
perc_limits = 0.05; % extra margin on the Y axis in percent with respect to the maximum margin of the datas
alpha_max=0.4;      % Increase the max Y margin so the plot is under the legend

          
%% Plot inizialization
fig = figure('name','Experimental results','Position',[Position, Dimension]);
set(fig,'defaulttextinterpreter','latex');  % latex compiler
set(fig,'DefaultTextFontname', 'CMU Time'); % font style

ha = tight_subplot(dim_raw,dim_col,marg_axes,marg_hBox,marg_wBox); 
x_limits = [data_3.time_window(1); data_3.time_window(end)];    % X limits

i = 1;

%% ========================================================================
%% p0 - X
axes(ha(i))
    set(gca,'FontSize',FontSize_axes) 
    plot(data_1.time_window, data_1.p0(:,1),'Color', r, 'LineWidth', LineWidth);
    hold on;
    plot(data_2.time_window, data_2.p0(:,1),'Color', g, 'LineWidth', LineWidth);
    plot(data_3.time_window, data_3.p0(:,1),'Color', b, 'LineWidth', LineWidth);
    line_fewer_markers(data_3.time_window, data_3.p0_d(:,1),10,'--*','Color', p, 'LineWidth', LineWidth, 'markersize', MarkerSize);
    set(gca,'Xtick',x_limits(1):xTick_step:x_limits(end))
    xlim(x_limits);
    grid on; 
    hold off;
    set(gca,'FontSize',FontSize_axes)        
    ylabel('$[\rm m]$','Interpreter','LaTex','FontSize', FontSize_label);
    datas_packed = padcat(data_1.p0(:,1), data_2.p0(:,1), data_3.p0(:,1), data_3.p0_d(:,1));
    y_limits = ylimits(datas_packed, perc_limits,alpha_max);
    ylim(y_limits);
    LegLocation = 'NorthEast'; %locationLegend(y_limits, eta(end,3)*r2g);
    h_leg = legend('$x_0^1$','$x_0^2$','$x_0^3$','$x_0^d$');
    set(h_leg,'orientation','horizontal','Location',LegLocation);
    set(h_leg,'FontSize',FontSize,'Interpreter','LaTex');     
    set(h_leg,'position',get(h_leg,'position') + [0 0.01 0 0]); 
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);

    index = get_figAxis(dim_fig,i);    
    ylabh(index(1),index(2)) = get(gca,'YLabel');
    ylabelPos = get(ylabh(index(1),index(2)),'Position');
    ylabelPosX(index(1),index(2),:) = ylabelPos;
    i = i+1;  
%% p0 - Z
axes(ha(i))
    k = 3;
    set(gca,'FontSize',FontSize_axes) 
    plot(data_1.time_window, data_1.p0(:,3),'Color', r, 'LineWidth', LineWidth);
    hold on;
    plot(data_2.time_window, data_2.p0(:,3),'Color', g, 'LineWidth', LineWidth);
    plot(data_3.time_window, data_3.p0(:,3),'Color', b, 'LineWidth', LineWidth);
    line_fewer_markers(data_3.time_window, data_3.p0_d(:,3),10,'--*','Color', p, 'LineWidth', LineWidth, 'markersize', MarkerSize);
    set(gca,'Xtick',x_limits(1):xTick_step:x_limits(end))
    xlim(x_limits);
    grid on; 
    hold off;
    set(gca,'FontSize',FontSize_axes)        
    ylabel('$[\rm m]$','Interpreter','LaTex','FontSize', FontSize_label);
    datas_packed = padcat(data_1.p0(:,k), data_2.p0(:,k), data_3.p0(:,k), data_3.p0_d(:,k));
    y_limits = ylimits(datas_packed, perc_limits,alpha_max);
    ylim(y_limits);
    LegLocation = 'NorthEast'; %locationLegend(y_limits, eta(end,3)*r2g);
    h_leg = legend('$z_0^1$','$z_0^2$','$z_0^3$','$z_0^d$');
    set(h_leg,'orientation','horizontal','Location',LegLocation);
    set(h_leg,'FontSize',FontSize,'Interpreter','LaTex');
    set(h_leg,'position',get(h_leg,'position') + [0 0.01 0 0]);      
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    
    index = get_figAxis(dim_fig,i);    
    ylabh(index(1),index(2)) = get(gca,'YLabel');
    ylabelPos = get(ylabh(index(1),index(2)),'Position');
    ylabelPosX(index(1),index(2),:) = ylabelPos;
    i = i+1;
    
%% ========================================================================
%% EE - x
axes(ha(i))
    k = 1;
    set(gca,'FontSize',FontSize_axes) 
    plot(data_1.time_window_arm, data_1.p_EE(:,k),'Color', r, 'LineWidth', LineWidth);
    hold on;
    plot(data_2.time_window_arm, data_2.p_EE(:,k),'Color', g, 'LineWidth', LineWidth);
    plot(data_3.time_window_arm, data_3.p_EE(:,k),'Color', b, 'LineWidth', LineWidth);
    line_fewer_markers(data_3.time_window, data_3.p_EE_d(:,k),10,'--*','Color', p, 'LineWidth', LineWidth, 'markersize', MarkerSize);
    set(gca,'Xtick',x_limits(1):xTick_step:x_limits(end))
    xlim(x_limits);
    grid on; 
    hold off;
    set(gca,'FontSize',FontSize_axes)        
    ylabel('$[\rm m]$','Interpreter','LaTex','FontSize', FontSize_label);
    datas_packed = padcat(data_1.p_EE(:,k), data_2.p_EE(:,k), data_3.p_EE(:,k), data_3.p_EE_d(:,k));
    y_limits = ylimits(datas_packed, perc_limits,alpha_max);
    ylim(y_limits);
    LegLocation = 'NorthEast'; %locationLegend(y_limits, eta(end,3)*r2g);
    h_leg = legend('$x_{ee}^1$','$x_{ee}^2$','$x_{ee}^3$','$x_{ee}^d$');
    set(h_leg,'orientation','horizontal','Location',LegLocation);
    set(h_leg,'FontSize',FontSize,'Interpreter','LaTex');  
    set(h_leg,'position',get(h_leg,'position') + [0 0.01 0 0]);
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    
    index = get_figAxis(dim_fig,i);    
    ylabh(index(1),index(2)) = get(gca,'YLabel');
    ylabelPos = get(ylabh(index(1),index(2)),'Position');
    ylabelPosX(index(1),index(2),:) = ylabelPos;
    i = i+1;
    
%% EE - z
axes(ha(i))
    k = 3;
    set(gca,'FontSize',FontSize_axes) 
    plot(data_1.time_window_arm, data_1.p_EE(:,k),'Color', r, 'LineWidth', LineWidth);
    hold on;
    plot(data_2.time_window_arm, data_2.p_EE(:,k),'Color', g, 'LineWidth', LineWidth);
    plot(data_3.time_window_arm, data_3.p_EE(:,k),'Color', b, 'LineWidth', LineWidth);
    line_fewer_markers(data_3.time_window, data_3.p_EE_d(:,k),10,'--*','Color', p, 'LineWidth', LineWidth, 'markersize', MarkerSize);
    set(gca,'Xtick',x_limits(1):xTick_step:x_limits(end))
    xlim(x_limits);
    grid on; 
    hold off;
    set(gca,'FontSize',FontSize_axes)        
    ylabel('$[\rm m]$','Interpreter','LaTex','FontSize', FontSize_label);
    datas_packed = padcat(data_1.p_EE(:,k), data_2.p_EE(:,k), data_3.p_EE(:,k), data_3.p_EE_d(:,k));
    y_limits = ylimits(datas_packed, perc_limits,alpha_max);
    ylim(y_limits);
    LegLocation = 'NorthEast'; %locationLegend(y_limits, eta(end,3)*r2g);
    h_leg = legend('$z_{ee}^1$','$z_{ee}^2$','$z_{ee}^3$','$z_{ee}^d$');
    set(h_leg,'orientation','horizontal','Location',LegLocation);
    set(h_leg,'FontSize',FontSize,'Interpreter','LaTex');  
    set(h_leg,'position',get(h_leg,'position') + [0 0.01 0 0]);    
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    
    index = get_figAxis(dim_fig,i);    
    ylabh(index(1),index(2)) = get(gca,'YLabel');
    ylabelPos = get(ylabh(index(1),index(2)),'Position');
    ylabelPosX(index(1),index(2),:) = ylabelPos;
    i = i+1; 
%% ========================================================================
% Attitude - PITCH
axes(ha(i))
    k = 2;
    rpy_off = 3; %[deg]
    set(gca,'FontSize',FontSize_axes) 
    plot(data_1.time_window, data_1.rpy0(:,k)*r2d-rpy_off,'Color', r, 'LineWidth', LineWidth);
    hold on;
    plot(data_2.time_window, data_2.rpy0(:,k)*r2d-rpy_off,'Color', g, 'LineWidth', LineWidth);
    plot(data_3.time_window, data_3.rpy0(:,k)*r2d-rpy_off,'Color', b, 'LineWidth', LineWidth);
    line_fewer_markers(data_3.time_window, data_3.rpy0_d(:,k)*r2d,10,'--*','Color', p, 'LineWidth', LineWidth, 'markersize', MarkerSize);
    set(gca,'Xtick',x_limits(1):xTick_step:x_limits(end))
    xlim(x_limits);
    grid on; 
    hold off;
    set(gca,'FontSize',FontSize_axes)        
    ylabel('$[\rm ^\circ]$','Interpreter','LaTex','FontSize', FontSize_label);
    datas_packed = padcat(data_1.rpy0(:,k), data_2.rpy0(:,k), data_3.rpy0(:,k))*r2d-rpy_off;
    y_limits = ylimits(datas_packed, perc_limits,alpha_max);
    ylim(y_limits);
    LegLocation = 'NorthEast'; %locationLegend(y_limits, eta(end,3)*r2g);
    h_leg = legend('$\theta_0^1$','$\theta_0^2$','$\theta_0^3$','$\theta_0^d$');
    set(h_leg,'orientation','horizontal','Location',LegLocation);
    set(h_leg,'FontSize',FontSize,'Interpreter','LaTex');  
    set(h_leg,'position',get(h_leg,'position') + [0 0.01 0 0]);    
    xlabel('Time $[\rm s]$','Interpreter','LaTex','FontSize', FontSize_label);
    
    index = get_figAxis(dim_fig,i);    
    ylabh(index(1),index(2)) = get(gca,'YLabel');
    ylabelPos = get(ylabh(index(1),index(2)),'Position');
    ylabelPosX(index(1),index(2),:) = ylabelPos;
    i = i+1;

%% ========================================================================    
%% Torque - y
axes(ha(i))
    k = 2;
    set(gca,'FontSize',FontSize_axes) 
    plot(data_1.time_window, data_1.tau0_c(:,k),'Color', r, 'LineWidth', LineWidth);
    hold on;
    plot(data_2.time_window, data_2.tau0_c(:,k),'Color', g, 'LineWidth', LineWidth);
    plot(data_3.time_window, data_3.tau0_c(:,k),'Color', b, 'LineWidth', LineWidth);
    line_fewer_markers(data_3.time_window, data_3.tau0_d(:,k),10,'--*','Color', p, 'LineWidth', LineWidth, 'markersize', MarkerSize);
    set(gca,'Xtick',x_limits(1):xTick_step:x_limits(end))
    xlim(x_limits);
    grid on; 
    hold off;
    set(gca,'FontSize',FontSize_axes)        
    ylabel('$[\rm Nm]$','Interpreter','LaTex','FontSize', FontSize_label);
    datas_packed = padcat(data_1.tau0_c(:,k), data_2.tau0_c(:,k), data_3.tau0_c(:,k), data_3.tau0_d(:,k));
    y_limits = ylimits(datas_packed, perc_limits,alpha_max);
    ylim(y_limits);
    LegLocation = 'NorthEast'; %locationLegend(y_limits, eta(end,3)*r2g);
    h_leg = legend('$u_{r_y}^1$','$u_{r_y}^2$','$u_{r_y}^3$','$u_{r_y}^d$');
    set(h_leg,'orientation','horizontal','Location',LegLocation);
    set(h_leg,'FontSize',FontSize,'Interpreter','LaTex');      
    set(h_leg,'position',get(h_leg,'position') + [0 0.01 0 0]);
    xlabel('Time $[\rm s]$','Interpreter','LaTex','FontSize', FontSize_label);
    
    index = get_figAxis(dim_fig,i);    
    ylabh(index(1),index(2)) = get(gca,'YLabel');
    ylabelPos = get(ylabh(index(1),index(2)),'Position');
    ylabelPosX(index(1),index(2),:) = ylabelPos;
    i = i+1; 

%% Correct Ylabel position
for i = 1:dim_col
    [minXi, minXiIndex]  = min(ylabelPosX(:,i,1));
     for j=1:dim_raw
         set(ylabh(j,i),'Position',[minXi ylabelPosX(j,i,2) ylabelPosX(j,i,3)]);
     end
end

% 
%% Save
name_exp = 'print_test';
directory = './../images/';
if save_ck == 1
    set(gcf, 'PaperSize', [15 15]);
    set(gcf, 'PaperPositionMode', 'auto');
    print(fig,'-dpdf',strcat(directory,name_exp));
    system(['pdfcrop ',directory,name_exp,'.pdf ',directory,name_exp,'.pdf']);
end
%% ========================================================================
