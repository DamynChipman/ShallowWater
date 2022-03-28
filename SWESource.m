function source = SWESource(x, y, t, q)
    % source = SWESource(x, y, t, q) Computes to source term at [x,y] and
    % time t with state q
    global g;
    
    % Unpack q
    h = q(1);
    u = q(2) / q(1);
    v = q(3) / q(1);
    
%     source = zeros(3,1);
    
%     if (x < 0.1 && x > -0.1)
%         source(2) = -g*h*10;
%     else
%         source(2) = 0;
%     end
    
    source = [0; 0; 0];

end