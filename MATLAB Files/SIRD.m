function x=SIRD(N,R,demands,alphas, mu,v,tau,iters,gamma,epsilon,cong_func,delta_func,proj_func)
    x=zeros(N*(R+1),iters);

    stops=zeros(N,1);
    stop=false;
    for i=1:N
        m=(i-1)*(R+1);
        x(m+1)=demands(i);
    end
    k=2;
    while ~stop && k<=iters
        chi=cong_func(x(:,k-1),N,R);
        delta_J_i=delta_func(x(:,k-1),alphas,mu,v,tau,chi,R,N);
        x(:,k)=x(:,k-1)-gamma*delta_J_i;
        x(:,k)=proj_func(x(:,k),-1*gamma*delta_J_i,N,R,demands);

        for i=1:N
            m=(i-1)*(R+1);
            if norm(x(m+1:m+R+1,k)-x(m+1:m+R+1,k-1))<epsilon
                stops(i)=1;
            end
        end

        if all(stops)
            stop=true;
        end

        if stop==true
            break
        end

        k=k+1;

        if k>iters
            break
        end
    end

    if stop==true
        disp("success, completed in "+int2str(k-1)+" iterations")
        for t=k+1:iters
            x(:,t)=x(:,k);
        end
    else
        disp("No success")
    end

end