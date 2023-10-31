% Estimating y'' = cos(x)y' - sin(x)y, 0 <= x <= 3pi/2 with y(0)=1 and y(3pi/2)= 1/e

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

error = NaN;
for n = [10 20 40 80]
    [x, y] = ODE_solver(p, q, r, a, b, cL, cR, n);
    hold on;
    plot(x, y, 'o');
    y_ex = f(x);
    h = (b-a)/n;

    new_error = norm(y_ex(:)-y(:), 'inf');
    fprintf("n: %i, Error norm: %f, Error ratio: %f\n", n, new_error, error/new_error);
    error = new_error;
end

waitfor(pl);