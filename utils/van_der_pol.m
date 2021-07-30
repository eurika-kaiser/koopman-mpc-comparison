function f = van_der_pol( t,x,u )

mu = 2;
f = [  x(2,:);
       -x(1,:)-mu.*(x(1,:).^2-1).*x(2,:)+u];

end