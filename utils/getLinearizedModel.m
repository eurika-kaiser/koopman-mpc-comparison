function sys = getLinearizedModel(f_ud,x_loc,u_loc)

n = length(x_loc);
x = sym('x',[n 1]); syms u;
f_ud_sym = f_ud(0,x,u);

Jx = jacobian(f_ud_sym,x);
Ju = jacobian(f_ud_sym,u);

sys.A = double(subs(Jx,[x;u],[x_loc;u_loc])); 
sys.B = double(subs(Ju,[x;u],[x_loc;u_loc]));
sys.c = double(subs(f_ud_sym,[x;u],[x_loc;u_loc])) - sys.A*x_loc - sys.B*u_loc;