% uprbs = (2*myprbs(Nsim,0.5) - 1);
uprbs = (((myprbs(Nsim, 0.5) - 0) * (exp_params.umax - exp_params.umin)) / (1 - 0)) + exp_params.umin;
u_dt = @(i)(  uprbs(i+1) );
f_cont_d = @(t,xx)( f_ud(t,xx,u_dt(t)) );

x = x0;

% Delayed initial condition (assume random control input in the past)
xstart = [Cy*x ; NaN(data_options.nD*(ny+m),1)];
if data_options.nD > 0
    for i = 1:data_options.nD
        urand = (((rand(m, 1) - 0) * (exp_params.umax - exp_params.umin)) / (1 - 0)) + exp_params.umin; % 2*rand(m,1) - 1;
        xp = f_ud((data_options.nD-i+1)*deltaT,x,urand);
        xstart = [Cy*xp ; urand; xstart(1:end-ny-m)];
        x = xp;
    end
else
    urand = 2*rand(m,1) - 1;
    xp = f_ud(0,x,urand);
    xstart = [Cy*xp];
    x = xp;
end

% Local linearization
x0 = xp;
xloc = xp;
sysLoc = getLinearizedModel(f_ud,xloc,urand); % ignores explicit time dependency for now!

% Inital conditions
x_true = xp;
for iM = 1:Nmodels
    if strcmp(koopman_options(iM).model, 'DMDc') || strcmp(koopman_options(iM).model, 'DMDc-2')
        xKM{iM} = sysKM(iM).liftFun(xp);
    else
        xKM{iM} = sysKM(iM).liftFun(xstart(1:koopman_options(iM).n_zeta));
    end
end

% Simulation
for i = 0:Nsim-1
    
    % True dynamics
    x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i)) ];
    
    % Koopman predictors
    for iM = 1:Nmodels
        xKM{iM} = [xKM{iM} sysKM(iM).A*xKM{iM}(:,end) + sysKM(iM).B*u_dt(i)];
    end
    
    % Local linearization predictor
    xloc = [xloc sysLoc.A*xloc(:,end) + sysLoc.B*u_dt(i) + sysLoc.c];
end