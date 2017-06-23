% Least mean squares (LMS) algorithm
function [e, w] = LMS(d, u, mu, M)

%% Initialization
N = length(u);
w = zeros(M,1);
e = zeros(N,1);
%% The LMS loop
for n=M:N
    % data pair
    u_temp = u(n:-1:n-M+1);    
    % update the filter
    y = w'*u_temp;
    e(n) = d(n) - y;
    w = w + mu*u_temp*conj(e(n));
end
w=conj(w);
end

