function [x, y] = ODE_solver(p, q, r, a, b, cL, cR, n, w=1.5)
    %%
    %
    % Generate a tridiagonal matrix from the given functions to solve the
    % equation
    %
    % y'' = p(x)y'+ q(x)y + r(s), a <= x <= b, f(a) = cL, f(b) = cR
    %
    % This will create a system of equations of the form
    % 
    % A * y = rhs
    %
    % Input:
    %   p  - function handle for p(x)
    %   q  - function handle for q(x)
    %   r  - function handle for r(x)
    %   a  - left boundary
    %   b  - right boundary
    %   cL - left boundary condition
    %   cR - right boundary condition
    %   n  - number of subintervals
    %   w  - relaxation parameter for SOR
    %
    % Output:
    %   
    %  x - vector of x values
    %  y - solution vector
    %
    %%

    h = (b-a)/n;

    x = a:h:b;

    right  =-1 + h/2 * p(x(2:end-2));
    middle = 2 + h^2 * q(x(2:end-1));
    left   =-1 - h/2 * p(x(3:end-1));

    A = diag(right, 1) + diag(middle) + diag(left, -1);

    rhs = (-h^2 * r(x(2:end-1)))';
    rhs(1) = rhs(1) + (1 + h*p(x(2))/2) * cL;
    rhs(end) = rhs(end) + (1 - h*p(x(end-1))/2) * cR ;

    y = SOR(A, rhs, ones(n-1, 1), 1e-11, w);
    y = [cL; y; cR];
end