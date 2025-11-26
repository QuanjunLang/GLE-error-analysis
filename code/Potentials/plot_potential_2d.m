function plot_potential_2d(potentialInfo, x_range, y_range, resolution)
%PLOT_POTENTIAL_2D Visualizes the 2D potential function over a specified domain
%
% Inputs:
%   potentialInfo - struct with .u (scalar multiplier) and .G (potential gradient)
%   x_range       - [xmin, xmax]
%   y_range       - [ymin, ymax]
%   resolution    - number of grid points along each axis

    if nargin < 4
        resolution = 100;
    end

    [X, Y] = meshgrid(linspace(x_range(1), x_range(2), resolution), ...
                      linspace(y_range(1), y_range(2), resolution));

    Z = zeros(size(X));

    % Reconstruct potential from gradient numerically (since G = ∇U is given)
    % Use: U(x) ≈ 0.5 * (x - x0)' * R * (x - x0) = 0.5 * G(x)
    for i = 1:resolution
        for j = 1:resolution
            xvec = [X(i,j); Y(i,j)];
            grad = potentialInfo.U(xvec);       % G(x) = (x - x0) * R * (x - x0)'
            Z(i,j) = 0.5 * grad;                 % Assume symmetric quadratic form
        end
    end

    Z = potentialInfo.u * Z;  % Scale if needed
    
    
    % Plot
    figure;
    subplot(121)
    surf(X, Y, Z, 'EdgeColor', 'none');
    colormap parula;
    colorbar;
    title('2D Potential Surface');
    xlabel('x'); ylabel('y'); zlabel('U(x, y)');
    view(30, 40); % nice 3D angle

    subplot(122);hold on;
    % Plot contou
    contourf(X, Y, Z, 30, 'LineColor', 'none'); % 30 contour levels
    colormap parula;
    colorbar;
    axis equal tight;
    title('2D Potential Contour');
    xlabel('x'); ylabel('y');

    plot(potentialInfo.U_argmin(1), potentialInfo.U_argmin(2), 'r.')
end