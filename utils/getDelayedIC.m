function [x0,zeta0] = getDelayedIC(x0,f_ud,m,Cy,deltaT,data_options,exp_params)

ny = size(Cy,1);

% Initial condition for the delay-embedded state assuming random control
% inputs in the past
x = x0;
zeta0 = [Cy*x ; NaN(data_options.nD*(ny+m),1)];
for i = 1:data_options.nD
    upast = (((rand(m, 1) - 0) * (exp_params.umax - exp_params.umin)) / (1 - 0)) + exp_params.umin; 
    xp = f_ud((data_options.nD-i+1)*deltaT,x,upast);
    zeta0 = [Cy*xp ; upast ; zeta0(1:end-ny-m)];
    x = xp;
end
x0 = x;
