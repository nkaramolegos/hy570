function [xi, w] = RLS(d, u, delta, lambda, M)
%rls  Recursive Least Squares

% initalizations
w=zeros(M,1);
P=eye(M)/delta;

% u and d should be column vectors as usuals

% input signal length
N=length(u);

% initialize error vector
xi=zeros(N,1);

% RLS loops
for n=M:N
    % data pair
    u_temp = u(n:-1:n-M+1);
    v=P*u_temp;    
    k=(1/(lambda+u_temp'*v))*v;
    xi(n)=d(n)-w'*u_temp;
    w=w+k*conj(xi(n));
    P=(1/lambda)*(P-k*v');
end
w=conj(w);
end

