function J=total_cost(x,alphas,mu,v,tau,chi,R,N)
    x_0=zeros(N,1);
    for i=1:N
        m=(i-1)*(R+1);
        x_0(i)=x(m+1);
    end
    diag_mu=zeros(R,R);
    for r=1:R
        diag_mu(r,r)=mu(r);
    end
    diag_mu_stack=diag_mu;
    chi_stack=chi;
    tau_stack=tau(1)*ones(R,1);
    v_stack=v;
    for i=2:N
        m=(i-1)*(R);
        diag_mu_stack=[diag_mu_stack,zeros(m,R);zeros(R,m),diag_mu];
        chi_stack=[chi_stack;chi];
        v_stack=[v_stack;v];
        tau_stack_i=tau(i)*ones(R,1);
        tau_stack=[tau_stack;tau_stack_i]; 
    end
    J=alphas'*x_0 + (diag_mu_stack*chi_stack+v_stack-tau_stack)'*(diag_mu_stack*chi_stack+v_stack-tau_stack);

end
