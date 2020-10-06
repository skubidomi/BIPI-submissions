function [gb_r, gb_i] = get_gabor_kernel(theta, sigma, gamma, lambda)
    [x, y] = meshgrid(-20:20);
    x_th = x*cos(theta)+y*sin(theta);
    y_th = -x*sin(theta)+y*cos(theta);
    
    gb_r = exp(-(x_th.^2 + gamma^2 * y_th.^2)/(2*sigma^2)) .* cos(2*pi*x_th/lambda);
    gb_i = exp(-(x_th.^2 + gamma^2 * y_th.^2)/(2*sigma^2)) .* sin(2*pi*x_th/lambda);
end

