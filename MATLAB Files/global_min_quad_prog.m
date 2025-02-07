function x=global_min_quad_prog(alphas,mu,v,tau,R,N,demands)
    diag_mu=zeros(R+1,R+1);
    for r=1:R
        diag_mu(r+1,r+1)=mu(r);
    end
    diag_mu_stack_1=diag_mu;
    for i=2:N
        diag_mu_stack_1=[diag_mu_stack_1,diag_mu]; 
    end
    diag_mu_stack=diag_mu_stack_1;
    alpha_stack=zeros(2*(R+1),1);
    alpha_stack(1)=alphas(1);
    v_stack=[0;v];
    tau_stack=[0;tau(1)*ones(R,1)];
    demand_stack=demands(1)*eye(R+1);
    for i=2:N
        m=(i-1)*(R+1);
        alpha_stack(m+1)=alphas(i);
        diag_mu_stack=[diag_mu_stack;diag_mu_stack_1];
        v_stack=[v_stack;0;v];
        tau_stack_i=[0;tau(i)*ones(R,1)];
        tau_stack=[tau_stack;tau_stack_i];
        demand_stack=[demand_stack,zeros(m,R+1);zeros(R+1,m),demands(i)*eye(R+1)];
    end
    % Define Q and c for the quadratic objective
    Q = 2 * (diag_mu_stack' * diag_mu_stack);
    c = alpha_stack + 2 * diag_mu_stack' * (v_stack - tau_stack);

    
    % Equality constraints
    A_eq = [ones(1,3), zeros(1,3); zeros(1,3), ones(1,3)];
    b_eq = demands;

    % Inequality constraints
    lb = zeros(6, 1); % x >= 0
    
    ub = demand_stack * ones(6,1); % upper bounds on x

    % Solve using quadprog
    options = optimoptions('quadprog', 'Display', 'iter');
    x = quadprog(Q, c, [], [], A_eq, b_eq, lb, ub, [], options);
    
    % Display solution
    % disp('Optimal solution:');
    % disp(x);


end

