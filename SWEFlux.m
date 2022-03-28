function flux = SWEFlux(q)
    %flux = SWEFlux(q) Computes the SW flux of q
    global g
    
    % Unpack q
    h = q(1);
    u = q(2) / q(1);
    v = q(3) / q(1);
    
    % Compute flux
    flux = [
        h*u, h*v;
        h*u*u + 0.5*g*h*h, h*u*v;
        h*u*v, h*v*v + 0.5*g*h*h;
    ];

end