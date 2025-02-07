function [x,NE_Check,no_NE,n_iters]=quadopti_BR_dynamics (alphas,mu,v,tau,R,N,demands,iters,epsilon,cong_hat_func, indiv_BR_func)
    x=zeros(N*(R+1),iters);

    stops=zeros(N,1);
    NE_Check=false;
    no_NE=false;
    for i=1:N
        m=(i-1)*(R+1);
        x(m+1)=demands(i);
    end
    diag_mu=zeros(R+1,R+1);
    for r=1:R
        diag_mu(r+1,r+1)=mu(r);
    end
    v_stack=[0;v];

    k=2;
    while ~NE_Check && k<=iters && ~no_NE
        for i=1:N
            m=(i-1)*(R+1);
            chi_hat_i=cong_hat_func(x(:,k-1),N,R,i);
            x(m+1:m+R+1,k)=indiv_BR_func(alphas,diag_mu,v_stack,tau,R,demands,chi_hat_i,i);
            if norm(x(m+1:m+R+1,k)-x(m+1:m+R+1,k-1))<epsilon
                stops(i)=1;
            end
        end

        if all(stops)
            NE_Check=true;
        end

        if NE_Check==true
            n_iters=k-1;
            break
        end

        if k>2
            if x(:,k)==x(:,k-2)
                no_NE=true;
            end
        end

        if no_NE==true
            n_iters=k-1;
            break
        end

        k=k+1;

        if k>iters
            n_iters=k-1;
            break
        end
    end

    if NE_Check==true
        disp("success, completed in "+int2str(k-1)+" iterations")
        for t=k+1:iters
            x(:,t)=x(:,k);
        end
    elseif no_NE==true
        disp("No success, solution has no NE")
        for t=k+1:iters
            x(:,t)=x(:,k);
        end
    else
        disp("No Success")
    end
end