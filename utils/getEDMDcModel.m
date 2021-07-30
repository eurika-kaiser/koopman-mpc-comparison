function [sys] = getEDMDcModel(X,Y,U,koopman_options)

%% Basis functions
basisFunction = koopman_options.basisFunction;
Nrbf = koopman_options.Nrbf;
n_zeta = koopman_options.n_zeta;

cent = rand(n_zeta,Nrbf)*2 - 1; % RBF centers
rbf_type = 'thinplate';
theta_max = pi;
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
Nlift = Nrbf + n_zeta;


%% Lift
disp('Starting LIFTING')

Xlift = liftFun(X);
Ylift = liftFun(Y);


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
sys.Cd = zeros(1,Nlift); sys.Cd(1) = 1;
sys.liftFun = liftFun;
sys.Nlift = Nlift;
fprintf('Regression for A, B, C DONE \n');

% Residual
fprintf( 'Regression residual %f \n', norm(Ylift - Alift*Xlift - Blift*U,'fro') / norm(Ylift,'fro') );
