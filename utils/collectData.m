function [X,Y,U,n_zeta,z_limits] = collectData(f_ud,Cy,data_options)
disp('Starting data collection')

n = data_options.n;
m = data_options.m;
Ntraj = data_options.Ntraj; % Number of trajectories
Nsim = data_options.Nsim; % Number of timesteps in each trajectory
nD = data_options.nD; % Number of delays
ny = size(Cy,1); % Number of outputs

xmin = data_options.xmin';
xmax = data_options.xmax';
umin = data_options.umin;
umax = data_options.umax;

% Random control input forcing
Ubig = (((rand(Nsim, Ntraj) - 0) * (data_options.umax - data_options.umin)) / (1 - 0)) + data_options.umin;
if isfield(data_options,'forcing')
    if strcmp(data_options.forcing,'sphs')
        for i = 1:Ntraj
            forcing = @(x,t) (2*randn.*sin(randn*t).*cos(t/(randn*10))).^2;
            u = forcing(0,[0:data_options.deltaT:(Nsim-1)*data_options.deltaT]);
            u = (((u - min(u)) * (data_options.umax - data_options.umin)) / (max(u) - min(u))) + data_options.umin;
            Ubig(:,i) = u;
        end
    end
end

% Random initial condition
if ~isfield(data_options,'X0')
    Xcurrent = (rand(n,Ntraj)*2 - 1);
else
    Xcurrent = data_options.X0;
end
X = []; Y = []; U = [];
zeta_current = [Cy*Xcurrent ; NaN(nD*(ny+m),Ntraj)];

n_zeta = (nD+1)*ny + nD*m; % dimension of the delay-embedded "state"
for i = 1:Nsim
    Xnext = f_ud(i*data_options.deltaT,Xcurrent,Ubig(i,:));
    
    % Clean-up & re-initialize 
    if any(isnan(Xnext(:)))
        idx = [];
        for iV = 1:n
            idx = [idx find(isnan(Xnext(iV,:)))];
        end
        idx = unique(idx);
        Xnext(:,idx) = (rand(n,length(idx))*2 - 1);
    end
    
    if isfield(data_options,'ylimits')
        if any(Xnext(:)<data_options.ylimits(1)) || any(Xnext(:)>data_options.ylimits(2))
            idx = [];
            for iV = 1:n
                idx = [idx find( Xnext(iV,:)<data_options.ylimits(1))];
                idx = [idx find( Xnext(iV,:)>data_options.ylimits(2))];
            end
            idx = unique(idx);
            Xnext(:,idx) = (rand(n,length(idx))*2 - 1);
        end
    end
    
    DX = Xnext - Xcurrent;
    if isfield(data_options,'dylimit')
        if any(abs(DX(:))>data_options.dylimit) 
            idx = [];
            for iV = 1:n
                idx = [idx find( DX(iV,:)<data_options.dylimit)];
            end
            idx = unique(idx);
            Xnext(:,idx) = (rand(n,length(idx))*2 - 1);
        end
    end
        
    % Update
    zeta_prev = zeta_current;
    if nD > 0
        zeta_current = [[Cy*Xnext ; Ubig(i,:)] ; zeta_current( 1:end-ny-m , : ) ];
    else
        zeta_current = Cy*Xnext;
    end
    if(i > nD) || nD==0 
        X = [X zeta_prev];
        Y = [Y zeta_current]; % yk, [uk-1, yk-1, uk-2, yk-2]
        U = [U Ubig(i,:)];
    end
    Xcurrent = Xnext;
    
end

% Create limit vector for scaling
z_limits = [Cy*xmin, Cy*xmax];
for i=1:nD
    z_limits = [z_limits; umin, umax; Cy*xmin, Cy*xmax];
end
    
fprintf('Data collection DONE \n');



function u = sphs(Pf,K,t)
% Schroeder-phased harmonic sequence
% Pf = 1; % Fundamental period
% K = 10;
u = zeros(size(t));
for i = 1:K
    theta = 2*pi/K * sum([1:i]);
    u = u + sqrt(2/(K))*cos(2*pi * i*t/Pf + theta) ;
end
