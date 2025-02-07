function delta_J_i=grad_indiv_cost(x,alphas,mu,v,tau,chi,R,N)
    alpha_stack=zeros(2*(R+1),1);
    alpha_stack(1)=alphas(1);
    diag_mu=zeros(R+1,R+1);
    for r=1:R
        diag_mu(r+1,r+1)=mu(r);
    end
    diag_mu_stack=diag_mu;
    chi_stack=[0;chi];
    v_stack=[0;v];
    tau_stack=[0;tau(1)*ones(R,1)];
    for i=2:N
        m=(i-1)*(R+1);
        alpha_stack(m+1)=alphas(i);
        diag_mu_stack=[diag_mu_stack,zeros(m,R+1);zeros(R+1,m),diag_mu];
        chi_stack=[chi_stack;0;chi];
        v_stack=[v_stack;0;v];
        tau_stack_i=[0;tau(i)*ones(R,1)];
        tau_stack=[tau_stack;tau_stack_i];
    end
    delta_J_i=alpha_stack+2*diag_mu_stack*(diag_mu_stack*chi_stack+v_stack-tau_stack);
end