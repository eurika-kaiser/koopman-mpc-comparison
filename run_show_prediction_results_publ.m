
if length(sysKM)>1
    xlimits = [min([min(Cy*x_true,[],2), min(Cy*xloc,[],2), min(sysKM(1).Cd*xKM{1},[],2), min(sysKM(2).Cd*xKM{2},[],2), min(sysKM(3).Cd*xKM{3},[],2)] ,[],2), ...
        max([max(Cy*x_true,[],2), max(Cy*xloc,[],2), max(sysKM(1).Cd*xKM{1},[],2), max(sysKM(2).Cd*xKM{2},[],2), max(sysKM(3).Cd*xKM{3},[],2)] ,[],2)];
else
    xlimits = [min([min(Cy*x_true,[],2), min(Cy*xloc,[],2), min(sysKM(1).Cd*xKM{1},[],2)],[],2), ...
        max([max(Cy*x_true,[],2), max(Cy*xloc,[],2), max(sysKM(1).Cd*xKM{1},[],2)],[],2),];
end
ulimits = [min(u_dt(0:Nsim-1)), max(u_dt(0:Nsim-1))];
du = 0.05*(ulimits(:,2)-ulimits(:,1));

if strcmp(select_system,'van_der_pol')
    cutoff_factor = 1.5;
else
    cutoff_factor = 2.5;
end
for i = 1:ny
    if xlimits(i,1)<cutoff_factor*min(Cy(i,:)*x_true)
        xlimits(i,1) = [cutoff_factor*min(Cy(i,:)*x_true)];
    end
    if xlimits(i,2) > cutoff_factor*max(Cy(i,:)*x_true)
        xlimits(i,2) = [cutoff_factor*max(Cy(i,:)*x_true)];
    end
end
dx = 0.05*(xlimits(:,2)-xlimits(:,1));
if ny>1
    dx(2) = 0.01*(xlimits(2,2)-xlimits(2,1));
end

ypos_x = -0.1;

%% Plot control input
figure,
subplot('Position',[0.15 0.25 0.8 0.7]); hold on, box on
stairs([0:Nsim-1]*deltaT,u_dt(0:Nsim-1),'k','linewidth',params_plotting.linewidth);
xlim([0 Tmax]), ylim([ulimits(1)-du, ulimits(2)+du])
set(gca,'FontSize',params_plotting.fontsize,'linewidth',params_plotting.linewidth_box);
xlbh = xlabel('Time','Interpreter','latex','Units','normalized');
ylbh = ylabel('Input','Interpreter','latex','Units','normalized');
ylbh.Position(1) = ypos_x;
set(gca,'TickLabelInterpreter','latex')

set(gcf,'Position',[100 100 600 280])
set(gcf,'PaperPositionMode','auto')
drawnow
pause(0.05)
print('-depsc2', '-loose', [pathfigs,['EX_',select_system,'_Prediction_Input','.eps']]);


%% Output
figure,
subplot('Position',[0.15 0.25 0.8 0.7]), box on
lw_koop = params_plotting.linewidth_models;
p1 = plot([0:Nsim]*deltaT,Cy*x_true,params_plotting.true.linestyle,'Color',params_plotting.true.color,'linewidth', lw_koop); hold on
for iM = 1:Nmodels
    ph{iM} = plot([0:Nsim]*deltaT,sysKM(iM).Cd*xKM{iM}, koopman_options(iM).linestyle, 'Color', koopman_options(iM).color,'linewidth',lw_koop);
end
p2 = plot([0:Nsim]*deltaT,Cy*xloc, params_plotting.loc.linestyle,'Color',params_plotting.loc.color,'linewidth',lw_koop-1);

xlabel('Time','Interpreter','latex','Units','normalized'),
ylbh = ylabel('Output','Interpreter','latex','Units','normalized'); ylbh.Position(1) = ypos_x;
axis([0 Tmax -4 4])
switch select_system
    case 'motor'
        axis([0 Tmax  xlimits])
    case 'van_der_pol'
        axis([0 Tmax min(xlimits(:,1)) max(xlimits(:,2))])
    case 'duffing'
        axis([0 Tmax xlimits])
end
set(gca,'FontSize',params_plotting.fontsize,'linewidth',params_plotting.linewidth_box);
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Position',[100 100 600 280])
set(gcf,'PaperPositionMode','auto')
drawnow
pause(0.05)

print('-depsc2', '-loose', [pathfigs,['EX_',select_system,'_Prediction_Output','.eps']]);
drawnow
pause(1)


%%
if exist('Error')~=0
    xtxt{1} = 'L-MPC';
    
    for iM = 1:Nmodels
        xtxt{1+iM} = koopman_options(iM).model;
    end
    
    MeanError = squeeze(mean(abs(Error),2));
    ErrorMedian = median( MeanError ,2);
    ErrorMin = min( MeanError ,[],2);
    ErrorMax = max( MeanError ,[],2);

    figure,
    subplot('Position',[0.15 0.35 0.8 0.6]), box on, hold on
    
    bhandle = bar(1, ErrorMedian(1));
    bhandle.FaceColor = params_plotting.loc.color;
    bhandle.EdgeColor = params_plotting.loc.color;
    for iM = 1:Nmodels
        bhandle = bar(iM+1, ErrorMedian(iM+1));
        bhandle.FaceColor = koopman_options(iM).color;
        bhandle.EdgeColor = koopman_options(iM).color;
    end
    er = errorbar(1:Nmodels+1,ErrorMedian,ErrorMin,ErrorMax);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    er.LineWidth = 3;
    set(gca,'yscale','log','FontSize',params_plotting.fontsize,'linewidth',params_plotting.linewidth_box);
    xlim([0.4 Nmodels+1+0.6])    
    ylbh = ylabel('Error','Interpreter','latex','Units','normalized'); ylbh.Position(1) = ypos_x;
    set(gca,'ytick',[10^0 10^1 10^2 10^3 10^4 10^5 10^6 10^7],'xtick',[1:Nmodels+1]+0.4,'xticklabels',xtxt,'TickLabelInterpreter','latex')
    xtickangle(25)
    set(gcf,'Position',[100 100 600 280])
    set(gcf,'PaperPositionMode','auto')
    drawnow
    pause(0.05)
    print('-depsc2', '-loose', [pathfigs,['EX_',select_system,'_Prediction_Stats','.eps']]);
end