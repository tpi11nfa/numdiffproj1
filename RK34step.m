%RK34step takes a single step with the classical RK4 method and puts 
%the result in unew. It uses the embedded RK3 to compute a local error 
%estimate which is stored in err.
function [unew, err] = RK34step(f, uold, told, h)
Y1 = f(told, uold);
Y2 = f(told + h/2, uold + h*Y1/2);
Y3 = f(told + h/2, uold + h*Y2/2);
Z3 = f(told + h, uold + h*(2*Y2 - Y1));
Y4 = f(told + h, uold + h*Y3);
unew = uold + h/6*(Y1 + 2*Y2 + 2*Y3 + Y4);
err = norm(h/6*(2*Y2 + Z3 - 2*Y3 - Y4));
end