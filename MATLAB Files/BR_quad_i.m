function x_i=BR_quad_i (alphas,diag_mu,v_stack,tau,R,demands,chi_hat_i,i)
    alpha_stack=[alphas(i);zeros(R,1)];
    tau_stack=[0;tau(i)*ones(R,1)];

    % Define Q and c for the quadratic objective
    Q = 2 * (diag_mu' * diag_mu);
    c = alpha_stack + 2 * diag_mu' * (diag_mu*[0;chi_hat_i] + v_stack - tau_stack);

    
    % Equality constraints
    A_eq = ones(1,R+1);
    b_eq = demands(i);

    % Inequality constraints
    lb = zeros(R+1, 1); % x >= 0
    
    ub = demands(i) * ones(R+1,1); % upper bounds on x

    % Solve using quadprog
    options = optimoptions('quadprog', 'Display', 'iter');
    x_i = quadprog(Q, c, [], [], A_eq, b_eq, lb, ub, [], options);
end