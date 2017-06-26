function [ e, w] = NLMS(d, u, beta, M, epsilon)
%NLMS Normalized LMS, like most of papers

%% Initialization
N = length(u);
w = zeros(M,1);
e = zeros(N,1);

%% The NLMS loop
for n=M:N
    % data pair
    u_temp = u(n:-1:n-M+1);    
    % update the filter
    y = w'*u_temp;
    e(n) = d(n) - y;
    w = w + (beta*u_temp)/(norm(u_temp).^2+epsilon)*conj(e(n));
end
w=conj(w);

end

