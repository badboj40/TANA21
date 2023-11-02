% Estimating y''=cos(x)y'-sin(x)y, 0<=x<=3pi/2 with y(0)=1 and y(3pi/2)=1/e

% Parameters
p = @(x) cos(x);
q = @(x) -sin(x);
r = @(x) zeros(size(x));
a = 0;
b = 3*pi/2;
cL = 1;
cR = 1/exp(1);

f = @(x) exp(sin(x));

x_reference = linspace(a, b, 1000);
y_reference = f(x_reference);
pl = plot(x_reference, y_reference, 'b');

n = [10 20 40 80];
w = [1.57 1.75 1.86 1.93];
error = NaN;
for i = 1:4
    [x, y] = ODE_solver(p, q, r, a, b, cL, cR, n(i), w(i));
    hold on;
    plot(x, y, 'diamond');
    y_ex = f(x);
    h = (b-a)/n(i);

    new_error = norm(y_ex(:)-y(:), 'inf');
    fprintf("n: %i, Error: %f, Ratio: %f\n", n(i), new_error, error/new_error);
    error = new_error;
end

waitfor(pl);