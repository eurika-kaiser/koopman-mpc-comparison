function J = eval_quadratic_cost_function(y, u, Q, R, QN, yr, ulin, qlin)

Nsim = size(u,2);
J = zeros(1,Nsim);
ny = size(y,1);
nu = size(u,1);

if nargin == 6
    ulin = zeros(nu,1); qlin = zeros(ny,1);
end

current_cost = @(yr,y,u) [ (y - yr)'*Q*(y - yr) + u'*R*u + ulin'*u + qlin'*y ];

for i = 1:Nsim-1
    J(i) = current_cost(yr(:,i),y(:,i),u(:,i));
end
J(Nsim) = (y(:,Nsim) - yr(:,Nsim))'*QN*(y(:,Nsim) - yr(:,Nsim));


