function run_example(select_system, pathfigs, show_figs_publ)

%% PLOTTING PARAMETERS
params_plotting.constraints.color = 'k';
params_plotting.constraints.linestyle = '-';
params_plotting.ref.color = [0.8,0.1,0.4];
params_plotting.ref.linestyle = '-';
params_plotting.true.color = 0.7*ones(1,3);
params_plotting.true.linestyle = '-';
params_plotting.loc.color = [0.2,0.8,0.2];
params_plotting.loc.linestyle = '--';
params_plotting.edmdc.color = 'b';
params_plotting.edmdc.linestyle = '-';
params_plotting.delaydmdc.color = [1,0.7,0.1]; 
params_plotting.delaydmdc.linestyle = '--';
params_plotting.delaydmdcvar.color = [1,0.4,0.1]; 
params_plotting.delaydmdcvar.linestyle = '-.';
params_plotting.dmdc.color = [0.6,0.1,0.7]; 
params_plotting.dmdc.linestyle = ':';
params_plotting.dmdcvar.color = [0.8,0.6,0.7];
params_plotting.dmdcvar.linestyle = '-.';
params_plotting.linewidth_models = 6;
params_plotting.linewidth = 4;
params_plotting.linewidth_box = 2;
params_plotting.fontsize = 28;



%% DYNAMICS

% PARAMETERS FOR DIFFERENT SYSTEMS 

switch select_system
    case 'motor'
        randgen_seed = 115123;
        f_u = @dyn_motor_scaled;
        n = 2; % Number of states
        m = 1; % Number of control inputs
        Cy = [0 1]; % Output matrix: y = Cy*x
        deltaT = 0.01;
        exp_params = struct('system', 'motor', ...
            'xmin',[-1,-1], 'xmax', [1,1], ...
            'Tmax_pred',1, 'Tmax_ctrl', 3, 'Tmax', 1, ...
            'REF', 'cos',...
            'umin', -1, 'umax', 1);
        data_options.nD = 16;
        
        
    case 'van_der_pol' % Van der Pol
        randgen_seed = 259278;
        f_u = @van_der_pol;
        n = 2; % Number of states
        m = 1; % Number of control inputs
        Cy = [0 1]; % Output matrix: y = Cy*x
        deltaT = 0.1;
        exp_params = struct('system', 'van_der_pol', ...
            'xmin',[-2,-2], 'xmax', [2,2], ... 
            'Tmax_pred',5, 'Tmax_ctrl', 10, 'Tmax', 3, ...
            'REF', 'sin', ...
            'umin', -5, 'umax', 5);
        data_options.nD = 8; 
        
        
    case 'duffing' % Duffing
        randgen_seed = 285123;
        f_u = @duffing;
        n = 2; % Number of states
        m = 1; % Number of control inputs
        Cy = [0 1]; % Output matrix: y = Cy*x
        deltaT = 0.1;
        exp_params = struct('system', 'duffing', ...
            'xmin',[-1.5,-1], 'xmax', [1.5,1], ...
            'Tmax_pred',5, 'Tmax_ctrl', 10, 'Tmax', 4, ...
            'REF', 'sin', ...
            'umin', -0.5, 'umax', 0.5);
        data_options.nD = 8; 
        
end

% 4th order Runge-Kutta 
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

k1bw = @(t,x,u) ( -f_u(t,x,u) );
k2bw = @(t,x,u) ( -f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3bw = @(t,x,u) ( -f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4bw = @(t,x,u) ( -f_u(t,x + k1(t,x,u)*deltaT,u) );
f_bw = @(t,x,u) ( x + (deltaT/6) * ( k1bw(t,x,u) + 2*k2bw(t,x,u) + 2*k3bw(t,x,u) + k4bw(t,x,u)  )   );

%% Collect data
rng(randgen_seed)
ny = size(Cy,1); % Number of outputs

data_options.Ntraj = 200;
data_options.Nsim = 1000;
data_options.n = n;
data_options.m = m;
data_options.deltaT = deltaT;
data_options.umin = exp_params.umin;
data_options.umax = exp_params.umax;
data_options.xmin = exp_params.xmin;
data_options.xmax = exp_params.xmax;

data_options.X0 = (( (rand(data_options.n, data_options.Ntraj) - 0) .* (exp_params.xmax - exp_params.xmin)') / (1-0)) + exp_params.xmin';
[X,Y,U,n_zeta,z_limits] = collectData(f_ud,Cy,data_options);

%% Koopman model regression
iModel = 1;
koopman_options(iModel).model = 'EDMDc';
koopman_options(iModel).nD = 1; 
if koopman_options(iModel).nD>data_options.nD
    koopman_options(iModel).nD=data_options.nD;
end
koopman_options(iModel).xmin = exp_params.xmin;
koopman_options(iModel).xmax = exp_params.xmax;
koopman_options(iModel).Nrbf = 100;
koopman_options(iModel).basisFunction = 'rbf';
koopman_options(iModel).rbf_type = 'thinplate'; 
koopman_options(iModel).output = [];
if ny>1
    koopman_options(iModel).output = [1:n];
end
koopman_options(iModel).n_zeta = koopman_options(iModel).nD*(ny+m)+ny;
koopman_options(iModel).color = params_plotting.edmdc.color;
koopman_options(iModel).linestyle = params_plotting.edmdc.linestyle;
koopman_options(iModel).z_limits = z_limits(1:koopman_options(iModel).n_zeta,:);
[sysKM(iModel)] = getKoopmanModel(X(1:koopman_options(iModel).n_zeta,:),Y(1:koopman_options(iModel).n_zeta,:),U,koopman_options(iModel)); %Alift, Blift, Clift, liftFun, koopman_options.Nlift


%% Linear Time delay model
iModel = iModel + 1;
koopman_options(iModel).basisFunction = 'tdc';
koopman_options(iModel).nD = 1; 
if koopman_options(iModel).nD>data_options.nD
    koopman_options(iModel).nD=data_options.nD;
end
koopman_options(iModel).model = ['DDMDc-',num2str(koopman_options(iModel).nD)];
koopman_options(iModel).xmin = exp_params.xmin;
koopman_options(iModel).xmax = exp_params.xmax;
koopman_options(iModel).n_zeta = koopman_options(iModel).nD*(ny+m)+ny;
koopman_options(iModel).color = params_plotting.delaydmdc.color;
koopman_options(iModel).linestyle = params_plotting.delaydmdc.linestyle;
koopman_options(iModel).output = [];
if ny>1
    koopman_options(iModel).output = [1:n];
end
koopman_options(iModel).z_limits = z_limits(1:koopman_options(iModel).n_zeta,:);
[sysKM(iModel)] = getKoopmanModel(X(1:koopman_options(iModel).n_zeta,:),Y(1:koopman_options(iModel).n_zeta,:),U,koopman_options(iModel));


if strcmp(select_system,'motor') || strcmp(select_system,'van_der_pol')
    iModel = iModel + 1;
    koopman_options(iModel).basisFunction = 'tdc';
    switch select_system
        case 'motor'
            koopman_options(iModel).nD = 16; 
        case 'van_der_pol'
            koopman_options(iModel).nD = 8; 
        case 'duffing'
            koopman_options(iModel).nD = 2;  
    end
    if koopman_options(iModel).nD>data_options.nD
        koopman_options(iModel).nD=data_options.nD;
    end
    koopman_options(iModel).model = ['DDMDc-',num2str(koopman_options(iModel).nD)];
    koopman_options(iModel).xmin = exp_params.xmin;
    koopman_options(iModel).xmax = exp_params.xmax;
    koopman_options(iModel).n_zeta = koopman_options(iModel).nD*(ny+m)+ny;
    koopman_options(iModel).color = params_plotting.delaydmdcvar.color;
    koopman_options(iModel).linestyle = params_plotting.delaydmdcvar.linestyle;
    koopman_options(iModel).output = [];
    if ny>1
        koopman_options(iModel).output = [1:n];
    end
    koopman_options(iModel).z_limits = z_limits(1:koopman_options(iModel).n_zeta,:);
    [sysKM(iModel)] = getKoopmanModel(X(1:koopman_options(iModel).n_zeta,:),Y(1:koopman_options(iModel).n_zeta,:),U,koopman_options(iModel));
end

%% Linear DMDc model
rng(randgen_seed)
data_options_dmdc = data_options;
data_options_dmdc.nD = 0;
[X,Y,U,~,z_limits] = collectData(f_ud,eye(n),data_options_dmdc);

iModel = iModel + 1;
koopman_options(iModel).model = 'DMDc';
koopman_options(iModel).basisFunction = 'linear';
koopman_options(iModel).n_zeta = n;
koopman_options(iModel).color = params_plotting.dmdc.color;
koopman_options(iModel).linestyle = params_plotting.dmdc.linestyle;
koopman_options(iModel).output = [];
koopman_options(iModel).xmin = exp_params.xmin;
koopman_options(iModel).xmax = exp_params.xmax;
koopman_options(iModel).z_limits = z_limits(1:koopman_options(iModel).n_zeta,:);
[sysKM(iModel)] = getKoopmanModel(X(1:koopman_options(iModel).n_zeta,:),Y(1:koopman_options(iModel).n_zeta,:),U,koopman_options(iModel)); %Alift, Blift, Clift, liftFun, koopman_options.Nlift
sysKM(iModel).C = Cy;
sysKM(iModel).Cd = Cy*eye(ny);

%% COMPARISON OF PREDICTORS
Nmodels = length(koopman_options);
Tmax = exp_params.Tmax_pred;
Nsim = Tmax/deltaT;
Nruns = 300;
Nruns_pre = Nruns+300;
Error = zeros(Nmodels+1,Nsim,Nruns);
Xbw = zeros(n,data_options.nD+1,Nruns_pre);
X0runs = zeros(n,Nruns);

% Initial conditions
X0pre = (( (rand(n, Nruns_pre) - 0) .* (exp_params.xmax - exp_params.xmin)') ./ (1-0)) + exp_params.xmin';

% Backward integrate so that initial conditions are uniformly distributed in domain
for iRun = 1:Nruns_pre
    % Delayed initial condition (assume random control input in the past)
    x = X0pre(:,iRun);
    Xbw(:,1,iRun) = x;
    xstart = [Cy*x ; NaN(data_options.nD*(ny+m),1)];
    if data_options.nD > 0
        for i = 1:data_options.nD
            urand = (((rand(m, 1) - 0) * (exp_params.umax - exp_params.umin)) / (1 - 0)) + exp_params.umin;
            xp = f_bw((data_options.nD-i+1)*deltaT,x,urand);
            xstart = [Cy*xp ; urand; xstart(1:end-ny-m)];
            x = xp;
            Xbw(:,i+1,iRun) = x;
        end
    end
    X0bw(:,iRun) = x;
end
% Filter unconverged training data (inf/nans)
for i = 1:n
    remove_nans = isnan(X0bw(i,:));
    X0bw(:,remove_nans) = [];
    remove_infs = isinf(X0bw(i,:));
    X0bw(:,remove_infs) = [];
    X0bw(:,abs(X0bw(i,:))>1e2) = [];
end
X0 = X0bw;

% Compute prediction error for ensemble of initial states
iRun = 1;
keep_ic = [];
for counter = 1:Nruns_pre
    x0 = X0(:,counter);
    run_prediction
    while any(isnan(x0)) || any(isinf(abs(x0))) || any(abs(x0)>1e2)
        x0 = X0(:,counter);
        run_prediction
        counter = counter + 1;
    end
    
    X0runs(:,iRun) = x0;
    Error(1,:,iRun) = (Cy*xloc(:,2:end)-Cy*x_true(:,2:end))./(Cy*x_true(:,2:end));
    for iM = 1:Nmodels
        Error(iM+1,:,iRun) = (sysKM(iM).Cd*xKM{iM}(:,2:end)-Cy*x_true(:,2:end))./(Cy*x_true(:,2:end));
    end
    keep_ic = [keep_ic, counter];
    iRun = iRun + 1;

    if iRun == Nruns+1
        break;
    end
end

%% PREDICTION

% Initial conditions
switch select_system
    case 'motor'
        x0 = [-0.0455; 0.3580];
    case 'duffing'
        x0 = [1, 0.5]';
    case 'van_der_pol'
        x0 = [0.5, 0.03]';
end
run_prediction

%% PLOT PREDICTION RESULTS
if show_figs_publ
    run_show_prediction_results_publ
end


save(fullfile(pathfigs,[select_system,'_prediction.mat']));

%% CONTROL

Tmax = exp_params.Tmax_ctrl;
Nsim = Tmax/deltaT;

switch exp_params.system
    case 'motor'
        switch exp_params.REF
            case 'cos'
                Q = eye(ny);
                R = 0.1; 
                ymin = -0.4;
                ymax = 0.4;              
                
                % Reference 
                yref = @(i,Nsim) [0.5*cos(2*pi*[i] / Nsim)]; 
                
                % Prediction horizon & weights
                Tpred = 10*deltaT;
                for iM =1:Nmodels
                    mpcKM(iM).Np = round(Tpred / deltaT);
                    mpcKM(iM).Q = Q;
                    mpcKM(iM).R = R;
                end
                Np_loc = round(Tpred / deltaT); 
                Np_dmdc = round(Tpred / deltaT); 
        end
    case 'van_der_pol'
        switch exp_params.REF       
            case 'sin'
                Q = eye(ny);
                R = 0.01;
                ymin = -0.8;
                ymax = 0.8;
                
                omega = 2*pi;
                yref = @(i,Nsim) [-1.0*sin(omega/(2*pi)*([1:Nsim]*deltaT))];
                
                % Prediction horizon
                Tpred = 1;
                for iM = 1:Nmodels
                    mpcKM(iM).Np = round(Tpred / deltaT);
                    mpcKM(iM).Q = Q;
                    mpcKM(iM).R = R;
                end
                Np_loc = round(Tpred / deltaT);
                Np_dmdc = round(Tpred / deltaT);      
                
                
        end
    case 'duffing'
        switch exp_params.REF
            case 'sin'
                Q = eye(ny);
                R = 0.1;
                ymin = -0.8;
                ymax = 0.8;
                
                omega = 2*pi;
                yref = @(i,Nsim) [-1.0*sin(omega/(2*pi)*([1:Nsim]*deltaT))]; 
                
                % Prediction horizon
                Tpred = 1;
                for iM = 1:Nmodels
                   mpcKM(iM).Np =  round(Tpred / deltaT);
                   mpcKM(iM).Q = Q;
                   mpcKM(iM).R = R;
                end
                Np_loc =  round(Tpred / deltaT);
                Np_dmdc = round(Tpred / deltaT);
        end
        
end

% Setting up controller/solver
umin = exp_params.umin;
umax = exp_params.umax;
for iM = 1:Nmodels
    % Constraints
    lb_y = [ymin ; nan(sysKM(iM).Nlift-1,1)];
    ub_y = [ymax ; nan(sysKM(iM).Nlift-1,1)];
    switch koopman_options(iM).model
        case 'EDMDc'
            lb_y = [ymin ; nan(sysKM(iM).Nlift-1,1)];
            ub_y = [ymax ; nan(sysKM(iM).Nlift-1,1)];
            
            
        case 'TDC'
            lb_y = [ymin ; nan(sysKM(iM).Nlift-1,1)];
            ub_y = [ymax ; nan(sysKM(iM).Nlift-1,1)];
    end
    mpcKM(iM).lb_y = lb_y;
    mpcKM(iM).ub_y = ub_y;
    mpcKM(iM).lb_u = umin;
    mpcKM(iM).ub_u = umax;
    
    % Build Koopman MPC controller
    mpcKM(iM).solvempc  = getMPC(sysKM(iM).A,sysKM(iM).B,sysKM(iM).Cd,0,mpcKM(iM).Q,mpcKM(iM).R,mpcKM(iM).Q, ...
        mpcKM(iM).Np,mpcKM(iM).lb_u, mpcKM(iM).ub_u, ...
        mpcKM(iM).lb_y, mpcKM(iM).ub_y,'qpoases');
    
end

% Build linear MPC object for DMDc
for iM = 1:Nmodels
    if strcmp(koopman_options(iM).model,'DMDc')
        ss_dmdc = ss(sysKM(iM).A,sysKM(iM).B,sysKM(iM).Cd,0,deltaT);
        MV = struct('Min',mpcKM(iM).lb_u,'Max',mpcKM(iM).ub_u);
        OV = struct('Min',mpcKM(iM).lb_y,'Max',mpcKM(iM).ub_y);
        WMO = struct('OutputVariables', Q, 'ManipulatedVariables', R);
        mpc_dmdc = mpc(ss_dmdc,ss_dmdc.Ts,Np_dmdc,Np_dmdc,WMO,MV,OV);
        xMPC = mpcstate(mpc_dmdc);
    end
end

% Initialization
X0runs = zeros(n,Nruns);
J_runs = zeros(Nmodels+1,Nruns);
Jy_runs = zeros(Nmodels+1,Nruns);

% Initial condition
switch select_system
    case 'motor'
        x0 = [-0.0455; 0.3580];
    case 'duffing'
        x0 = [1, 0.5]';
    case 'van_der_pol'
        x0 = [0.5, 0.03]';
end        
[x0,zeta0] = getDelayedIC(x0,f_ud,m,Cy,deltaT,data_options,exp_params);

% Reference
yrr = yref([1:Nsim],Nsim);

% Run control simulations
run_control

% Plot control results
XX_dmdc = [];
if show_figs_publ
    run_show_control_results_publ
end

save(fullfile(pathfigs,[select_system,'_control.mat']));

