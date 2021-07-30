
x_loc = x0;
x_unc = x0;
x_dmdc = x0;
xMPC.Plant = x_dmdc;

XX_koop = zeros(n,Nsim,Nmodels);
UU_koop = zeros(m,Nsim,Nmodels);
koopmod_was_infeas = zeros(Nmodels,1);
for iM = 1:Nmodels
    x_koop{iM} = x0;
    zeta{iM} = zeta0(1:koopman_options(iM).n_zeta); % Delay-embedded "state"
    XX_koop(:,1,iM) = x0;
    koopmod_was_infeas(iM) = 0;
end
koop_ind_infeas = cell(Nmodels,1);

XX_loc = x0; UU_loc = [];
XX_unc = x0;
XX_dmdc = x0; UU_dmdc = [];

u_loc = 0;

wasinfeas= 0;
ind_inf = [];

% Closed-loop simultion start


for i = 0:Nsim-1
    if(mod(i,10) == 0)
        fprintf('Closed-loop simulation: iterate %i out of %i \n', i, Nsim)
    end
    
    % Current value of the reference signal
    yr = yrr(:,i+1);
    
    % Koopman MPC

    for iM = 1:Nmodels
        if strcmp(koopman_options(iM).model, 'DMDc') || strcmp(koopman_options(iM).model, 'DMDc-2')
            xlift = sysKM(iM).liftFun(zeta{iM}); % Lift
        else
            xlift = sysKM(iM).liftFun(zeta{iM}); % Lift
        end
        [u_koop{iM}, J_koop, flag_koop] = mpcKM(iM).solvempc(xlift,yr); % Get control input
        x_koop{iM} = f_ud(i*deltaT,x_koop{iM},u_koop{iM}); % Update true state; starts with t=0
        if strcmp(koopman_options(iM).model, 'DMDc') || strcmp(koopman_options(iM).model, 'DMDc-2')
            zeta{iM} = [ x_koop{iM}];
        else
            zeta{iM} = [ Cy*x_koop{iM} ; u_koop{iM}; zeta{iM}(1:end-ny-m)]; % Update delay-embedded state
        end
        if flag_koop~=0
            koopmod_was_infeas(iM) = 1;
            koop_ind_infeas{iM} = [koop_ind_infeas{iM} i];
        end
    end
    
    % DMDc with MPC object
    [u_dmdc,~] = mpcmove(mpc_dmdc,xMPC,[],yr);
    x_dmdc = f_ud(i*deltaT,x_dmdc,u_dmdc); % Update true state
    xMPC.Plant = x_dmdc;
    
    % Local linearization MPC
    sysLoc = getLinearizedModel(f_ud,x_loc,u_loc); %TODO add time dependency
    ub_y_loc = nan(n,1);
    lb_y_loc = nan(n,1);
    for iV = 1:ny
        ub_y_loc(Cy(iV,:)==1) = ymax;
        lb_y_loc(Cy(iV,:)==1) = ymin;
    end
    Qloc = eye(ny);
    [U_loc,~,optval] = solveMPCprob(sysLoc.A,sysLoc.B,Cy,sysLoc.c,Qloc,R,Qloc,Np_loc,umin, umax,lb_y_loc,ub_y_loc,x_loc,yr); % Get control input
    u_loc = U_loc(1:m,1);
    if(optval == Inf) % Detect infeasibility
        ind_inf = [ind_inf i];
        wasinfeas = 1;
    end
    x_loc = f_ud(i*deltaT,x_loc,u_loc); % Update true state
    
    % Uncontrolled
    x_unc = f_ud(i*deltaT,x_unc,0); % Update true state
    
    % Store values
    for iM = 1:Nmodels
        XX_koop(:,i+2,iM) = x_koop{iM};
        UU_koop(:,i+1,iM) = u_koop{iM};
    end
    XX_loc = [XX_loc x_loc];
    UU_loc = [UU_loc u_loc];
    XX_unc = [XX_unc x_unc];
    XX_dmdc = [XX_dmdc x_dmdc];
    UU_dmdc = [UU_dmdc u_dmdc];
end

if(isempty(ind_inf))
    ind_inf = Nsim;
end

x_unc = x0;
XX_unc = x0;
for i = 0:10*Nsim-1
    % Uncontrolled
    x_unc = f_ud(i*deltaT,x_unc,0);
    XX_unc = [XX_unc x_unc];
end

% Cost
yrr_cap = yrr;

J_loc = eval_quadratic_cost_function(Cy*XX_loc, UU_loc, Q, R, Q, yrr_cap);
J_koop = zeros(Nmodels,Nsim);
for iM = 1:Nmodels
    J_koop(iM,:) = eval_quadratic_cost_function(Cy*XX_koop(:,:,iM), UU_koop(:,:,iM), Q, R, Q, yrr_cap);
end

yrr_cap = yrr;
if ~isnan(ymin) || ~isempty(ymin)
    yrr_cap(yrr_cap<ymin) = ymin;
    yrr_cap(yrr_cap>ymax) = ymax;
end

Jy_loc = eval_output_error(Cy*XX_loc, yrr_cap);
Jy_koop = zeros(Nmodels,Nsim);
for iM = 1:Nmodels
    Jy_koop(iM,:) = eval_output_error(Cy*XX_koop(:,:,iM), yrr_cap);
end
