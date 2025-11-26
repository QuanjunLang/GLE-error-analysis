function plot_trajectory_highdim(trajInfo, sysInfo)
%PLOT_TRAJECTORY_HIGHDIM Plots high-dimensional trajectory over time
%
% Inputs:
%   trajInfo: struct with fields
%       - x [d × N]: position trajectory
%       - v [d × N]: velocity trajectory
%   sysInfo: struct with field
%       - dt: time step size

    x = trajInfo.x;
    v = trajInfo.v;
    [d, N] = size(x);
    t = sysInfo.t_grid;

    % Plot position
    figure;
    subplot(2,1,1);
    plot(t, x');
    xlabel('Time t');
    ylabel('Position');
    title('High-Dimensional Position Trajectories');
    legend(arrayfun(@(i) sprintf('x_%d', i), 1:d, 'UniformOutput', false));
    grid on;

    % Plot velocity
    subplot(2,1,2);
    plot(t, v');
    xlabel('Time t');
    ylabel('Velocity');
    title('High-Dimensional Velocity Trajectories');
    legend(arrayfun(@(i) sprintf('v_%d', i), 1:d, 'UniformOutput', false));
    grid on;
end
