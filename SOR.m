function [x,k] = SOR(A, b, x0, tol, w)
%%
%  Gauss-Seidel iterative method to approximate the solution of a
%  linear system Ax=b up to a user defined tolerance
%
%  INPUT: 
%    A   - n by n square, non-singular matrix
%    b   - n by 1 right hand side vector
%    x0  - n by 1 vector containing that initial guess for the iteration
%    tol - user set tolerance for the stopping condition in the iteration
%    w   - value should be 1 <= w < 2 
%
%  OUTPUT:
%    x - n by 1 vector containing the iterative solution
%    k - number of iterations
%%
%  get the system size
    n = length(A);
%%
%  Gauss-Seidel iteration which overwrites the current approximate solution
%  with the new approximate solution (pseudocode available in the lecture
%  notes on page 46)
    x = x0;
    max_its = 22000; % Max iterations allowed

    for k = 1:max_its
        x_old = x;
        for i = 1:n
            sum = 0;
            for j = 1:n
                if j ~= i
                    sum = sum + A(i,j) * x(j);
                end
            end    
            x(i) = (1 - w) * x_old(i) + w * ((b(i) - sum) / A(i, i));    
        end
        r = b - A*x;
        if norm(r, 2) <= tol * norm(b, 2)
            break
        end
    end
%
end