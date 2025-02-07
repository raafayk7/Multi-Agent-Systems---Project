function x_proj=proj(x,x_diff,N,R,demands)
    x_proj=zeros(2*(R+1),1);
    for i=1:N
        m=(i-1)*(R+1);
        for r=0:R
            x_diff(m+r+1)=x_diff(m+r+1)+max(0,-1*x(m+r+1))-max(0,x(m+r+1)-demands(i));
            x_proj(m+r+1)=x(m+r+1)+max(0,-1*x(m+r+1))-max(0,x(m+r+1)-demands(i));
        end
    end
    inverse_eye=ones(R+1,R+1)-eye(R+1);
    inverse_eye_stack=inverse_eye;
    for i=2:N
        m=(i-1)*(R+1);
        inverse_eye_stack=[inverse_eye_stack,zeros(m,R+1);zeros(R+1,m),inverse_eye];
    end
    x_proj=x_proj-(1/R)*inverse_eye_stack*x_diff;

end