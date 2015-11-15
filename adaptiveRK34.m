%combines RK34step and newstep to create an adaptive ODE solver which keeps
%the error estimate equal to tol with an adaptive step size
function [t,y] = adaptiveRK34(f, y0, t0, tf, tol)
h = (abs(tf - t0)*tol.^(1/4))/(100*(1 + norm(f(t0, y0))));
err = tol;
t(1) = t0;
y{1} = y0;
i = 1;
k = 4;
while t(end) <= tf
    errold = err;
    [y{i+1}, err] = RK34step(f, y{i}, t(i), h);
    i = i + 1;
    h = newstep(tol, err, errold, h, k);
    t(i) = t(i - 1) + h(1); 
end
h = tf - t(i - 1);
t(i) = tf;
[y{i},~] = RK34step(f, y{i - 1}, t(i - 1), h);
end