
xlimits = [min([min(XX_unc,[],2), min(XX_loc,[],2), squeeze( min(XX_koop,[],2))],[],2), max([max(XX_unc,[],2), max(XX_loc,[],2), squeeze( max(XX_koop,[],2))],[],2)];
for i = 1:ny
    xlimits(Cy(i,:)==1,:) = [min([yrr(i,:),Cy(i,:)*xlimits(:,1)],[],2), max([yrr(i,:),Cy(i,:)*xlimits(:,2)], [], 2)];
end
ulimits = [min([umin,min(UU_loc,[],2), squeeze( min(UU_koop,[],2))'],[],2), max([umax,max(UU_loc,[],2), squeeze( max(UU_koop,[],2))'],[],2)];

switch select_system
    case 'motor'
        xlimits(2,1) = -0.5;
    case 'van_der_pol'
        xlimits(2,:) = [-1.1, 1.1];
        xlimits(2,:) = [-1.1,1.1];
end

legtxt{1} = 'Truth/Uncontrolled';
legtxt{2} = 'Constraints';
legtxt{3} = 'Reference';
legtxt{4} = 'L-MPC';

for iM = 1:Nmodels
    legtxt{4+iM} = koopman_options(iM).model;
end

%% Control signal
du = 0.05*(ulimits(:,2)-ulimits(:,1));
figure,
subplot('Position',[0.15 0.25 0.8 0.7]), hold on, box on
rectangle('Position',[0+0.01,umax,Tmax-0-0.02,du(1)-0.1*du(1)], 'FaceColor',0.9*ones(1,3), 'EdgeColor',0.9*ones(1,3), 'LineWidth',0.5)
rectangle('Position',[0+0.01,umin-du(1)+0.1*du(1),Tmax-0-0.02,du(1)-0.1*du(1)], 'FaceColor',0.9*ones(1,3), 'EdgeColor',0.9*ones(1,3), 'LineWidth',0.5)
p3 = plot([0:Nsim]*deltaT,umax*ones(Nsim+1,1),params_plotting.constraints.linestyle,'Color',params_plotting.constraints.color,'linewidth',params_plotting.linewidth); hold on
p1 = plot([0:Nsim]*deltaT,umin*ones(Nsim+1,1),params_plotting.constraints.linestyle,'Color',params_plotting.constraints.color,'linewidth',params_plotting.linewidth); hold on
for iM = 1:Nmodels
    ph{iM} = plot([0:Nsim-1]*deltaT,UU_koop(:,1:Nsim,iM),koopman_options(iM).linestyle,'Color',koopman_options(iM).color,'linewidth',params_plotting.linewidth_models); hold on
end
p2 = plot([0:Nsim-1]*deltaT,UU_loc(1:Nsim),params_plotting.loc.linestyle,'Color',params_plotting.loc.color,'linewidth',params_plotting.linewidth_models-1); hold on
if ~isempty(XX_dmdc)
    plot([0:Nsim-1]*deltaT,UU_dmdc(:,1:Nsim),params_plotting.loc.linestyle,'Color','m','linewidth',params_plotting.linewidth_models-1);
end
axis([0,Tmax, ulimits(1)-du(1), ulimits(2)+du(1)] )
xlbh = xlabel('Time','Interpreter','latex','Units','normalized');
ylbh = ylabel('Input','Interpreter','latex','Units','normalized'); ylbh.Position(1) = ypos_x;
set(gca,'FontSize',params_plotting.fontsize, 'LineWidth',params_plotting.linewidth_box);
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Position',[100 100 600 280])
set(gcf,'PaperPositionMode','auto')
drawnow
pause(0.05)
print('-depsc2', '-loose', [pathfigs,['EX_',select_system,'_Control_Input','.eps']]);

pause(1)
%% Output
figure,
subplot('Position',[0.15 0.25 0.8 0.7]), hold on, box on
if ~isnan(ymax)
    rectangle('Position',[0+0.01,ymax,Tmax-0-0.02,xlimits(2,2)-ymax-0.01], 'FaceColor',0.9*ones(1,3), 'EdgeColor',0.9*ones(1,3), 'LineWidth',0.5)
end
if ~isnan(ymin)
    rectangle('Position',[0+0.01,xlimits(2,1)+0.01,Tmax-0-0.02,(ymin-xlimits(2,1)-0.01)], 'FaceColor',0.9*ones(1,3), 'EdgeColor',0.9*ones(1,3), 'LineWidth',0.5)
end
p1=plot([0:Nsim]*deltaT,ymax*ones(Nsim+1,1),params_plotting.constraints.linestyle,'Color',params_plotting.constraints.color,'linewidth',params_plotting.linewidth); hold on
p2=plot([0:Nsim]*deltaT,ymin*ones(Nsim+1,1),params_plotting.constraints.linestyle,'Color',params_plotting.constraints.color,'linewidth',params_plotting.linewidth);
p3=plot([0:Nsim-1]*deltaT, yrr,params_plotting.ref.linestyle,'Color',params_plotting.ref.color,'linewidth',params_plotting.linewidth);
for iM = 1:Nmodels
    ph{iM} = plot([0:Nsim-1]*deltaT,Cy*XX_koop(:,1:Nsim,iM),koopman_options(iM).linestyle,'Color',koopman_options(iM).color,'linewidth',params_plotting.linewidth_models); hold on
end
p4 = plot([0:Nsim-1]*deltaT,Cy*XX_loc(:,1:Nsim),params_plotting.loc.linestyle,'Color',params_plotting.loc.color,'linewidth',params_plotting.linewidth_models-1);
if ~isempty(XX_dmdc)
    plot([0:Nsim-1]*deltaT,Cy*XX_dmdc(:,1:Nsim),params_plotting.loc.linestyle,'Color','m','linewidth',params_plotting.linewidth_models-1);
end
axis([0,Tmax,min(Cy*xlimits(:,1),[],1),max(Cy*xlimits(:,2),[],1)])
set(gca,'FontSize',params_plotting.fontsize, 'LineWidth',params_plotting.linewidth_box);
ylbh = ylabel('Output','Interpreter','latex','Units','normalized'); ylbh.Position(1) = ypos_x;
xlabel('Time','Interpreter','latex','Units','normalized'),
set(gca,'TickLabelInterpreter','latex')

switch select_system
    case 'motor'
        ylabel('Output','Interpreter','latex','Position',ylbh.Position),
        xlabel('Time','Interpreter','latex','Position',xlbh.Position),
end
set(gcf,'Position',[100 100 600 280])
set(gcf,'PaperPositionMode','auto')
drawnow
pause(0.05)
print('-depsc2', '-loose', [pathfigs,['EX_',select_system,'_Control_Output','.eps']]);

pause(1)
%% Cost
jlimits = [min( [cumsum(J_loc); cumsum(J_koop,2)], [], 2)  max( [cumsum(J_loc); cumsum(J_koop,2)], [], 2)];
jlimits = [min(jlimits(:,1)), max(jlimits(:,2))];
jlimits(1) = 0;
maxvals = max( [cumsum(J_loc); cumsum(J_koop,2)], [], 2);
[~,IX] = sort(maxvals,'descend');
if maxvals(IX(1))/maxvals(IX(2)) > 5
    jlimits(2) = maxvals(IX(2));
end

dj = 0.05*(jlimits(2)-jlimits(1));

figure,
subplot('Position',[0.15 0.25 0.8 0.7]), hold on, box on
for iM = 1:Nmodels
    plot([0:Nsim-1]*deltaT,cumsum(J_koop(iM,:)),koopman_options(iM).linestyle,'Color',koopman_options(iM).color,'linewidth',params_plotting.linewidth_models)
end
plot([0:Nsim-1]*deltaT,cumsum(J_loc),params_plotting.loc.linestyle,'Color',params_plotting.loc.color,'linewidth',params_plotting.linewidth_models-1)
if ~isempty(XX_dmdc)
    plot([0:Nsim-1]*deltaT,cumsum(J_dmdc),params_plotting.loc.linestyle,'Color','m','linewidth',params_plotting.linewidth_models-1);
end
% return
axis([0,Tmax,jlimits(1),jlimits(2)+dj])
ylbh = ylabel('Cost','Interpreter','latex','Units','normalized'); ylbh.Position(1) = ypos_x;
xlabel('Time','Interpreter','latex','Units','normalized'),
set(gca,'FontSize',params_plotting.fontsize, 'LineWidth',params_plotting.linewidth_box);
set(gca,'TickLabelInterpreter','latex')

set(gcf,'Position',[100 100 600 280])
set(gcf,'PaperPositionMode','auto')
drawnow
pause(0.05)
print('-depsc2', '-loose', [pathfigs,['EX_',select_system,'_Control_Cost','.eps']]);

pause(1)

%% Tracking Error
jlimits = [min( [cumsum(Jy_loc); cumsum(Jy_koop,2)], [], 2)  max( [cumsum(Jy_loc); cumsum(Jy_koop,2)], [], 2)];
jlimits = [min(jlimits(:,1)), max(jlimits(:,2))];
jlimits(1) = 0;
maxvals = max( [cumsum(Jy_loc); cumsum(Jy_koop,2)], [], 2);
[~,IX] = sort(maxvals,'descend');
if maxvals(IX(1))/maxvals(IX(2)) > 5
    jlimits(2) = maxvals(IX(2));
end

dj = 0.05*(jlimits(2)-jlimits(1));

figure,
subplot('Position',[0.15 0.25 0.8 0.7]), hold on, box on
for iM = 1:Nmodels
    plot([0:Nsim-1]*deltaT,cumsum(Jy_koop(iM,:)),koopman_options(iM).linestyle,'Color',koopman_options(iM).color,'linewidth',params_plotting.linewidth_models)
end
plot([0:Nsim-1]*deltaT,cumsum(Jy_loc),params_plotting.loc.linestyle,'Color',params_plotting.loc.color,'linewidth',params_plotting.linewidth_models-1)
axis([0,Tmax,jlimits(1),jlimits(2)+dj])
ylbh = ylabel('Error','Interpreter','latex','Units','normalized'); ylbh.Position(1) = ypos_x;
xlabel('Time','Interpreter','latex','Units','normalized'),
set(gca,'FontSize',params_plotting.fontsize, 'LineWidth',params_plotting.linewidth_box);
set(gca,'TickLabelInterpreter','latex')

set(gcf,'Position',[100 100 600 280])
set(gcf,'PaperPositionMode','auto')
drawnow
pause(0.05)
print('-depsc2', '-loose', [pathfigs,['EX_',select_system,'_Control_Error','.eps']]);

pause(1)
%% Phase plot
figure,
subplot('Position',[0.15 0.25 0.8 0.7]), hold on, box on
dx = 0.005*(xlimits(1,2)-xlimits(1,1));
dy = 0.01*(xlimits(2,2)-xlimits(2,1));

if ~isnan(ymax)
    rectangle('Position',[xlimits(1,1)+dx,ymax,xlimits(1,2)-xlimits(1,1)-2*dx,xlimits(2,2)-ymax-dy], 'FaceColor',0.9*ones(1,3), 'EdgeColor',0.9*ones(1,3), 'LineWidth',0.5)
end
if ~isnan(ymin)
    rectangle('Position',[xlimits(1,1)+dx,xlimits(2,1)+dy,xlimits(1,2)-xlimits(1,1)-2*dx,(ymin-xlimits(2,1)-dy)], 'FaceColor',0.9*ones(1,3), 'EdgeColor',0.9*ones(1,3), 'LineWidth',0.5)
end
if n==2
    p1 = plot(XX_unc(1,:),XX_unc(2,:),params_plotting.true.linestyle,'Color',params_plotting.true.color,'linewidth',params_plotting.linewidth_models-1); hold on
    p2 = plot([xlimits(1,1),xlimits(1,2)], [ymin,ymin],params_plotting.constraints.linestyle,'Color',params_plotting.constraints.color,'linewidth',params_plotting.linewidth); hold on
    p3 = plot([xlimits(1,1),xlimits(1,2)], [ymax,ymax],params_plotting.constraints.linestyle,'Color',params_plotting.constraints.color,'linewidth',params_plotting.linewidth); hold on
    if size(yrr,1)==2
        plot(yrr(1,:),yrr(2,:),params_plotting.ref.linestyle,'Color',params_plotting.ref.color,'linewidth',params_plotting.linewidth); hold on
    end
    for iM = 1:Nmodels
        ph{iM} = plot(XX_koop(1,:,iM),XX_koop(2,:,iM),koopman_options(iM).linestyle,'Color',koopman_options(iM).color,'linewidth',params_plotting.linewidth_models); hold on
    end
    p4 = plot(XX_loc(1,1:Nsim),XX_loc(2,1:Nsim),params_plotting.loc.linestyle,'Color',params_plotting.loc.color,'linewidth',params_plotting.linewidth_models-1);
    axis([xlimits(1,1), xlimits(1,2), xlimits(2,1), xlimits(2,2)])
elseif n==3
    p1 = plot3(XX_unc(1,:),XX_unc(2,:),XX_unc(3,:),params_plotting.true.linestyle,'Color',params_plotting.true.color,'linewidth',params_plotting.linewidth_models-1); hold on
    for iM = 1:Nmodels
        ph{iM} = plot3(XX_koop(1,:,iM),XX_koop(2,:,iM),XX_koop(3,:,iM),koopman_options(iM).linestyle,'Color',koopman_options(iM).color,'linewidth',params_plotting.linewidth_models); hold on
    end
    p4 = plot3(XX_loc(1,1:Nsim),XX_loc(2,1:Nsim),XX_loc(3,1:Nsim),params_plotting.loc.linestyle,'Color',params_plotting.loc.color,'linewidth',params_plotting.linewidth_models-1);
    
    if ~isempty(XX_dmdc)
        plot3(XX_dmdc(1,1:Nsim),XX_dmdc(2,1:Nsim),XX_dmdc(3,1:Nsim),params_plotting.loc.linestyle,'Color','m','linewidth',params_plotting.linewidth_models-1);
    end
    axis([xlimits(1,1), xlimits(1,2), xlimits(2,1), xlimits(2,2), xlimits(3,1), xlimits(3,2)])
    zlabel('$x_3$','Interpreter','latex'),
    view(-100,20) %view(25,20)
end
ylbh = ylabel('$x_2$','Interpreter','latex','Units','normalized'); ylbh.Position(1) = ypos_x;
xlabel('$x_1$','Interpreter','latex','Units','normalized'),
set(gca,'FontSize',params_plotting.fontsize, 'LineWidth',params_plotting.linewidth_box);
set(gca,'TickLabelInterpreter','latex')

set(gcf,'Position',[100 100 600 400])
set(gcf,'PaperPositionMode','auto')
drawnow
pause(0.05)
print('-depsc2', '-loose', [pathfigs,['EX_',select_system,'_Control_Phase','.eps']]);

pause(1)
%%
figure,
subplot('Position',[0.15 0.7 0.8 0.1]), hold on, box on
p1 = plot([0:Nsim]*deltaT,Cy*XX_unc(:,1:Nsim+1),params_plotting.true.linestyle,'Color',params_plotting.true.color,'linewidth',lw_koop); hold on
p2 = plot([0:Nsim]*deltaT,umax*ones(Nsim+1,1),params_plotting.constraints.linestyle,'Color',params_plotting.constraints.color,'linewidth',lw_koop-1); hold on
p3 = plot([0:Nsim-1]*deltaT, yrr,params_plotting.ref.linestyle,'Color',params_plotting.ref.color,'linewidth',lw_koop);
p4 = plot([0:Nsim]*deltaT,Cy*XX_loc(:,1:Nsim+1),params_plotting.loc.linestyle,'Color',params_plotting.loc.color,'linewidth',lw_koop);
for iM = 1:Nmodels
    ph{iM} = plot([0:Nsim]*deltaT,Cy*XX_koop(:,:,iM),koopman_options(iM).linestyle,'Color',koopman_options(iM).color,'linewidth',lw_koop); hold on
end
if ~isempty(XX_dmdc)
    plot([0:Nsim]*deltaT,Cy*XX_dmdc(:,:),params_plotting.loc.linestyle,'Color','m','linewidth',params_plotting.linewidth_models-1);
end
axis off

phandles = [p1(1),p2(1),p3(1),p4(1)];
for iM = 1:Nmodels
    phandles = [phandles,ph{iM}(1)];
end

LEG  = legend(phandles, legtxt,'box','off');
set(LEG,'Interpreter','latex','orientation','horizontal','position',[0.1 0.25 0.8 0.05],'fontsize',params_plotting.fontsize+2)
set(gcf,'Position',[100 100 1500 150])
set(gcf,'PaperPositionMode','auto')
drawnow
pause(0.05)
print('-depsc2', '-loose', [pathfigs,['EX_',select_system,'_Control_Legend','.eps']]);


