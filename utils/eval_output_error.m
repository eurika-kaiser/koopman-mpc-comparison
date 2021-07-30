function J = eval_output_error(y, yr)

Nsim = size(yr,2);
J = zeros(1,Nsim);

current_cost = @(yr,y) [ (y - yr)'*(y - yr) ];

for i = 1:Nsim
    J(i) = current_cost(yr(:,i),y(:,i));
end


