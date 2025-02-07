function chi_hat_i=cong_minus_i(x,N,R,i)
    chi_hat_i=zeros(R,1);
    for r=1:R
        for j=1:N
            if j~=i
                m=(j-1)*(R+1);
                chi_hat_i(r)=chi_hat_i(r)+x(m+r+1);
            end
        end
    end
end