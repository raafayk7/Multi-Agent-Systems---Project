function J_i=indiv_cost(x,alphas,mu,v,tau,chi,N,R)
    J_i=zeros(N,1);
    for i=1:N
        m=(i-1)*(R+1);
        tau_stack=tau(i)*ones(R,1);
        diag_mu=zeros(R,R);
        for r=1:R
            diag_mu(r,r)=mu(r);
        end
        J_i(i)=alphas(i)*x(m+1) + (diag_mu*chi + v - tau_stack)'*(diag_mu*chi + v - tau_stack);
    end
end