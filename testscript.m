%%
%task 1.1
%test RK4step for the linear test equation and plot the error in a log-log
%diagram. plot next to t,t^4 in linlog to verify that error is of order 4
%(i.e. that they have the same slope)
ltest = @(t,y) -0.2*y;
h = zeros(1, 7);
errors = zeros(1, 7);
for i = 2:8
    N = 2.^i;
    t = linspace(0, 1, N + 1);
    y = zeros(1, N + 1);
    err = zeros(1, N + 1);
    y(1) = 1;
    h(i - 1) = 1/N;
    for k = 2:N
        y(k) = RK4step(ltest, y(k -1), t(k), h(i - 1));
        err(k) = y(k) - expm(t(k)*(-0.2));
    end
    errors(i - 1) = norm(err);
end
loglog(h, errors, '*')
grid on
hold on
ycomp = h.^4;
plot(h, ycomp)
xlabel('log h','FontSize',12)
ylabel('log err','FontSize',12)
legend('= log error', ' = h^4')
hold off

%%
%task 1.4: test adaptive RK34 method
ltest = @(t,y) -0.2*y;
y0 = 1;
t0 = 0;
tf = 10;
tol = 1e-6;
[t, y] = adaptiveRK34(ltest, y0, t0, tf, tol);
plot(t,cell2mat(y), '*')
hold on
plot(t, exp(-0.2*t))
plot(t,0)

hold off

%%
%task 2: lotka-volterra equation, solution periodicity
tol = 1e-8;
A = [3 9 15 15];
y0 = [1; 1];
t0 = 0;
tf = 12;
%dudt represents the time derivative of the lotka-volterra equation system
dudt = @(t,u) [A(1)*u(1) - A(2)*u(1)*u(2); A(3)*u(1)*u(2) - A(4)*u(2)];
[t,y] = adaptiveRK34(dudt, y0, t0, tf, tol);
y = cell2mat(y);
plot(t,y)
xlabel('t','FontSize',12)
ylabel('x,y','FontSize',12)
figure
plot(y(1,:), y(2,:))
xlabel('x','FontSize',12)
ylabel('y','FontSize',12)

%%
%task 2: check i numerical H(x,y) stays near initial value
tol = 1e-6;
A = [3 9 15 15];
y0 = [1; 1];
t0 = 0;
tf = 100;
%dudt represents the time derivative of the lotka-volterra equation
dudt = @(t,u) [A(1)*u(1) - A(2)*u(1)*u(2); A(3)*u(1)*u(2) - A(4)*u(2)];
H = @(u) A(3)*u(1) + A(2)*u(2) - A(4)*log(u(1)) - A(1)*log(u(2));
[t,y] = adaptiveRK34(dudt, y0, t0, tf, tol);
H0 = H(y0);
Hdiff = @(u) abs(H(u)/H0 - 1);
ylen = length(y);
Hd = zeros(1, ylen);
for i = 1:ylen
    Hd(i) = Hdiff(y{i});
end
semilogy(t, Hd)
xlabel('t','FontSize',12)
ylabel('numerical vs analytical H','FontSize',12)

%%
%task 3.1
mu = 100;
%dudt is the time derivative of the van der pol equation
dudt = @(t,u) [u(2); mu*(1 - u(1).^2)*u(2) - u(1)];
tol = 1e-8;
y0 = [1; 0.5];
t0 = 0;
tf = 200;
[t,y] = adaptiveRK34(dudt, y0, t0, tf, tol);
y = cell2mat(y);
plot(t,y)
xlabel('t','FontSize',12)
ylabel('y1, y2','FontSize',12)
figure
plot(y(1,:), y(2,:))
xlabel('y1','FontSize',12)
ylabel('y2','FontSize',12)


%%
%task 3.2
mu = [10 15 22 33 47 68 100 150 220 330 470 680 1000];
%according to the plot, its stiffness changes between mu = 470 and 680,
%so at around mu = 500, it starts being stiff
mulen = length(mu);
ts = zeros(1, mulen);
tol = 1e-6;
y0 = [2; 0];
t0 = 0;
%dudt is the time derivative of the van der pol equation
for i = 1:mulen
    tf = 0.7*mu(i);
    dudt = @(t,u) [u(2); mu(i)*(1 - u(1).^2)*u(2) - u(1)];
    [t,y] = adaptiveRK34(dudt, y0, t0, tf, tol);
    ts(i) = length(t);
end
loglog(mu, ts, '*')
xlabel('log mu','FontSize',12)
ylabel('log N','FontSize',12)
hold on
semilogy(mu, mu.^2) %show that N ~ C*mu^2 (i.e. q = 2)
hold off

%%
%task 3.3
mu = 10;
dudt = @(t,u) [u(2); mu*(1 - u(1).^2)*u(2) - u(1)];
tol = 1e-4;
y0 = [2; 0];
t0 = 0;
tf = 200;
ts = [t0 tf];%TSPAN for ode15s
tic;
%[t,y] = adaptiveRK34(dudt, y0, t0, tf, tol);
'adaptiveRK34'
toc
tic;
[t1, y1] = ode15s(dudt, ts, y0);
'ode15s'
toc
y = cell2mat(y);
plot(t,y, '*')
figure
plot(t1,y1)
%ode15s gets slightly faster as mu gets larger, and outclasses adaptiveRK34
%for small mu also if tolerance is small. step size difference is very visible for smaller
%tolerance levels like 1e-6, and adaptiveRK34 doesn't get quicker for small
%mu until larger tolerance values like 1e-4

%%
mu = [10 15 22 33 47 68 100 150 220 330 470 680 1000 10000];
y0 = [2; 0];
t0 = 0;
tf = 200;
ts = [t0 tf];
mulen = length(mu);
ti = zeros(1, mulen);
for i = 1:mulen
        dudt = @(t,u) [u(2); mu(i)*(1 - u(1).^2)*u(2) - u(1)];
        [t1, y1] = ode15s(dudt, ts, y0);
        ti(i) = length(t1);
end
loglog(mu, ti, '*')
xlabel('log mu','FontSize',12)
ylabel('log N','FontSize',12)
        
        
        
