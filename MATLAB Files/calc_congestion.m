function chi=calc_congestion(x,N,R)
    chi=zeros(R,1);
    for r=1:R
        for i=1:N
            m=(i-1)*(R+1);
            chi(r)=chi(r)+x(m+r+1);
        end
    end
end