% Finding optimal w-values for N = [10, 20, 40, 80]
for N=[10 20 40 80]
    A = full(gallery('tridiag',N,-1,2,-1));
    known_x = ones(N,1); % make known solution vector of 1s
    manuf_b = A*known_x; % manufacture the right-hand-side
    x0 = zeros(N, 1);
    M = 100;
    w = zeros(M, 1);
    iterations = zeros(M, 1);
    
    for i = 1:M
        w_i = 1 + (i-1)/M;
        [x,k] = SOR(A, manuf_b, x0, 1e-11, w_i);
        w(i) = w_i;
        iterations(i) = k;
    end
    
    
    semilogy(w, iterations);
    hold on;
    legend('N = 10', 'N = 20', 'N = 40', 'N = 80')
    
    [minVal, k] = min(iterations);
    
    fprintf("N: %i, Optimal w: %f\n", N, w(k));
end