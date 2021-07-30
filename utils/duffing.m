function f = duffing( t,x,u )


a = 1; b = -1; d = -0.3; f0 = 0.5; omega = 1.2; 
f = [  x(2,:);
       a*x(1,:)+b*x(1,:).^3+d*x(2,:)+f0*cos(omega*t)+(3+cos(x(1,:))).*u];

end 
