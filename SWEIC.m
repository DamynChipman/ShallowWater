function q0 = SWEIC(x, y)
    % q0 = SWEIC(x, y) Computes the initial condition for the SWE
    hu0 = 0;
    hv0 = 0;
%     if (x > 0 && y > 0)
%         h0 = 1;
%     else
%         h0 = 0.5;
%     end

    gaussian = @(x, y, A, x0, y0, sigma_x, sigma_y) A*exp(-(((x - x0)^2)/(2*sigma_x^2) + ((y - y0)^2)/(2*sigma_y^2)));

%     h0 = gaussian(x, y, A, x0_ll, y0_ll, sigma_x, sigma_y) + ...
%          gaussian(x, y, A, x0_ur, y0_ur, sigma_x, sigma_y);
    
    n_peaks = 10;
    x0 = 0;
    ys = linspace(-0.5, 0.5, n_peaks);
    A = 6 / n_peaks;
    sigma_x = 0.3;
    sigma_y = 0.3;
    h0 = 0.4;
    for i=1:5
        h0 = h0 + gaussian(x, y, A, x0, ys(i), sigma_x, sigma_y);
    end

    % Random noise in IC
    noise_amplitude = 0.2;
    h0 = h0 + noise_amplitude * rand(1);

%     h0 = 0.5;
    q0 = [h0; hu0; hv0];
    
end