function [sys] = getKoopmanModel(X,Y,U,koopman_options)

%% Basis functions
basisFunction = koopman_options.basisFunction;
n_zeta = koopman_options.n_zeta;

switch basisFunction
    case 'linear'
        liftFun = @(xx)( [xx] );
        Nlift = n_zeta;
        
        sys = lift_and_regress(X,Y,U,liftFun);
        sys.Cd = eye(Nlift);
        sys.liftFun = liftFun;
        sys.Nlift = Nlift;
        
    case 'tdc'
        liftFun = @(xx)( [xx] );
        Nlift = n_zeta;
        
        sys = lift_and_regress(X,Y,U,liftFun);
        sys.Cd = zeros(1,Nlift); sys.Cd(1) = 1;
        sys.liftFun = liftFun;
        sys.Nlift = Nlift;
        
    case 'rbf'
        % Basis functions
        Nrbf = koopman_options.Nrbf;
        cent = (( (rand(n_zeta, Nrbf) - 0) .* (koopman_options.z_limits(:,2) - koopman_options.z_limits(:,1))) / (1-0)) + koopman_options.z_limits(:,1);
        rbf_type = koopman_options.rbf_type;
        theta_max = pi;
        liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
        Nlift = Nrbf + n_zeta;
        
        sys = lift_and_regress(X,Y,U,liftFun);
        sys.Cd = zeros(1,Nlift); sys.Cd(1) = 1;
        sys.liftFun = liftFun;
        sys.Nlift = Nlift;       
        
end
fprintf('Regression for A, B, C DONE \n');

if ~isempty(koopman_options.output)
    sys.Cd = zeros(length(koopman_options.output),Nlift); 
    for i = 1:length(koopman_options.output)
        sys.Cd(i,koopman_options.output(i)) = 1;    
    end
end

function sys = lift_and_regress(X,Y,U,liftFun)
%% Lift
disp('Starting LIFTING')

Xlift = liftFun(X);
Ylift = liftFun(Y);
Nlift = size(Xlift,1);

%% Regression

disp('Starting REGRESSION for A,B,C')

W = [Ylift ; X];
V = [Xlift ; U];
VVt = V*V';
WVt = W*V';
ABC = WVt * pinv(VVt);
Alift = ABC(1:Nlift,1:Nlift);
Blift = ABC(1:Nlift,Nlift+1:end);
if size(X,1)>1
    Clift = ABC(Nlift+1:end,1:Nlift);
else
    Clift = zeros(size(X,1),size(Alift,2));
    Clift(1:size(X,1),1:size(X,1)) = eye(size(X,1));
end

sys.A = Alift;
sys.B = Blift;
sys.C = Clift;

% Residual
fprintf( 'Regression residual %f \n', norm(Ylift - Alift*Xlift - Blift*U,'fro') / norm(Ylift,'fro') );