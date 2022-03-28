function lambda = SWELambda(q, n)
    % lambda = SWELambda(q, n) Computes the eigenvector of the SWE
    global g;

    % Unpack q
    h = q(1);
    u = q(2) / q(1);
    v = q(3) / q(1);
    
    % Compute wave speed
    c = sqrt(g * h);
    
    % Compute lambda
    un = u*n(1) + v*n(2);
    lambda = [un - c, un, un + c];
    
end