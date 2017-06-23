% Least Squares algorithm

function  h  = LS(u, d, L)

    % the matrix A for the trainnig sequence, constructed from the known to
    % receiver trainning sequence
    A=toeplitz(u,[u(1) zeros(1,L-1)]);  
    h=((A'*A)^-1)*A'*d;

end

