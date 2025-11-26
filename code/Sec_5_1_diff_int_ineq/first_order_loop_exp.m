close all
clear all
clc
addpaths

rng(1);

% dV_t = -gamma * V_t * dt - (integral from 0 to t of K(t, s) * V_s ds) * dt + sigma * dB_t

%% System settings
sysInfo.d       = 1;

sysInfo.dt      = 0.01;          % Time step
sysInfo.T       = 100;             % Total simulation time

sysInfo.sigma   = 0;
sysInfo.Batch   = 100;

sysInfo         = update_sysInfo(sysInfo);


all_gamma = linspace(1, 10, 10);
all_decay = linspace(1, 6, 12);

all_fitInfo = cell(length(all_gamma), length(all_decay));

for i = 1:length(all_gamma)
    for j = 1:length(all_decay)
        i 
        j
        gamma = all_gamma(i);
        decay = all_decay(j);
        amp1 = 1;

        sysInfo.gamma = all_gamma(i);

        kernel_type = 'exp';                       % Specify the type of kernel to generate ('osc', 'exp', 'powerlaw', 'diff')
        kernelInfo = kernel_info(kernel_type, sysInfo, 'decay', decay, 'amp1', amp1); % Generate the kernel functions
        % plot_kernels(kernelInfo);                       % Plot the kernel functions in three panels: linear, semilog-y, and log-log




        %% Generate trajectory
        noise = get_noise(sysInfo);

        traj_input{1}.v0 = randn(sysInfo.d, 1);
        traj_input{1}.noise = noise;

        trajInfo{1} = generate_data_first_order_highdim(sysInfo, kernelInfo.K1, traj_input{1});

        % traj_input{2}.v0 = randn(sysInfo.d, 1);
        % traj_input{2}.noise = noise;
        % trajInfo{2} = generate_data_first_order_highdim(sysInfo, kernelInfo.K2, traj_input{2});


        %%
        fitInfo = plot_velocity_highdim(trajInfo{1}.v, sysInfo, 'traj 1', false);
        % plot_velocity_highdim(trajInfo{2}.v, sysInfo, 'traj 2');
        % fitInfo = plot_velocity_highdim(abs(trajInfo{1}.v - trajInfo{2}.v), sysInfo, 'traj difference', false);

        all_fitInfo{i, j} = fitInfo;

        get_theoretical_rate_exponential(decay, gamma, amp1)
    end
end
%% collect fitted rate
for i = 1:length(all_gamma)
    for j = 1:length(all_decay)
        rate(i, j) = all_fitInfo{i, j}.a_semilogy;
    end
end


for i = 1:length(all_gamma)
    for j = 1:length(all_decay)
        
        gamma = all_gamma(i);
        beta = all_decay(j);
        
        rate_thm(i, j) = get_theoretical_rate_exponential(beta, gamma, amp1);
    end
end



%%

% Assuming:
% - all_gamma is a row or column vector of length G
% - all_decay is a vector of length D
% - rate is a G × D matrix where rate(i, j) corresponds to gamma(i), decay(j)

figure; hold on;

for j = 1:length(all_decay)
    plot(all_gamma, rate(:, j), 'DisplayName', sprintf('\\beta = %.2f', all_decay(j)));
end

xlabel('\gamma');
ylabel('Rate (a_{semilogy})');
legend('show');
title('Decay Rate vs \gamma for Different \beta');
grid on;


%%
% Assuming all_gamma and all_rate are vectors
[nGamma, nRate] = size(all_fitInfo);
rsquare_values = zeros(nGamma, nRate);

for i = 1:nGamma
    for j = 1:nRate
        rsquare_values(i, j) = all_fitInfo{i, j}.gof_semi.rsquare;
    end
end

% Plotting
figure;
subplot(131)
imagesc(all_gamma, all_decay, rsquare_values);  % Now x=gamma, y=decay
colorbar;
xlabel('\gamma (all\_gamma)');
ylabel('\beta (all\_rate)');
title('R-squared from Semilogy Fit');
set(gca, 'YDir', 'normal');  % Optional: adjust as needed

subplot(132)
imagesc(all_gamma, all_decay, rate);  % Now x=gamma, y=decay
colorbar;
xlabel('\gamma (all\_gamma)');
ylabel('\beta (all\_rate)');
title('Rate Fit');
set(gca, 'YDir', 'normal');  % Optional: adjust as needed


subplot(133)
imagesc(all_gamma, all_decay, rate_thm);  % Now x=gamma, y=decay
colorbar;
xlabel('\gamma (all\_gamma)');
ylabel('\beta (all\_rate)');
title('Rate Theory');
set(gca, 'YDir', 'normal');  % Optional: adjust as needed


%%
% 
% % Initialize figure
% figure;
% hold on;
% 
% % Plot rate vs gamma for each beta
% for i = 1:2:length(all_decay)
%     beta = all_decay(i);
% 
%     plot(all_gamma, rate_thm(:, i), 'LineWidth', 1.5, ...
%          'DisplayName', sprintf('\\beta = %.1f', beta));
% 
%     plot(all_gamma, rate(:, i), '-.', 'LineWidth', 1.5, ...
%          'DisplayName', sprintf('\\beta = %.1f', beta));
% end
% 
% xlabel('a');
% ylabel('\beta');
% % title('Theoretical decay rate and fitted decay rate');
% legend('Location', 'best');
% grid on;
% 
% 
% exportgraphics(gcf, '/Users/quanjunlang/Documents/GitHub/Tracjectory-prediction-error-of-GLE/fig/ExpDecayEx.pdf', 'ContentType','vector');



%%


% Choose ONE of these size presets (inches). Comment out the other.
figW = 3.8; figH = 3.8;   % two-column figure (papers, LaTeX)
% figW = 3.5; figH = 3.0; % single-column compact

fontSize = 12;            % 9–10 for papers, 12–14 for slides
lineW    = 1.0;           % overlay line width

% === Figure (tight layout) ===
fig = figure('Units','inches','Position',[0 0 figW figH], 'Color','w');
hold on;

% Colors
color_true = [0 0 0];    % black for theoretical
color_fit  = [0.85 0 0]; % dark red for fitted


plot_gamma = [1,2,3,4,5,6,8];


for k = 1:1:length(plot_gamma)
    i = plot_gamma(k);
    a = all_gamma(i);

    % Theoretical rate (thin solid, same color)
    h1 = plot(all_decay, rate_thm(i, :), 'LineWidth', 1.0, ...
        'LineStyle', '-', ...
        'Color', color_true);

    % Fitted rate (bold dashed, same color)
    h2 = plot(all_decay, rate(i, :), 'LineWidth', 2.0, ...
        'LineStyle', '--', ...
        'Color', color_fit);

    % Add text label near the end of the curve
    text(all_decay(end)+0.1, rate(i,end), sprintf('a=%d', a), ...
        'VerticalAlignment','middle','HorizontalAlignment','left', 'FontSize',fontSize);
end

plot(all_decay, -all_decay, 'LineWidth', 1.0, ...
        'LineStyle', '-', ...
        'Color', color_true);
text(all_decay(end)+0.1, -all_decay(end), 'a=\infty', ...
        'VerticalAlignment','middle','HorizontalAlignment','left', 'FontSize',fontSize,'Interpreter','tex');


xlabel('\beta'); 
xlim([1, 6])
title('Exponential decay rate $r$', 'Interpreter','latex');
grid on;

ax = gca;
ax.Position(3) = ax.Position(3) * 0.98;  % shrink width a bit to leave margin on right
set(gca, 'Position', get(gca, 'Position') + [0 0.02 0 -0.05])
% Single legend entry for line styles
legend([h1 h2], {'Theoretical','Fitted'}, 'Location','best');

set(gca, 'FontSize', fontSize)




exportgraphics(gcf, '/Users/quanjunlang/Documents/GitHub/Tracjectory-prediction-error-of-GLE/fig/ExpDecayEx.pdf', 'ContentType','vector');

%%
% % Define the range and resolution for beta and gamma
% beta_vals = linspace(0, 5, 100);    % adjust range as needed
% gamma_vals = linspace(0, 5, 100);   % adjust range as needed
% 
% % Create meshgrid
% [Beta, Gamma] = meshgrid(beta_vals, gamma_vals);
% 
% % Fixed constant
% C = 1;
% 
% % Evaluate the function over the grid
% Rate = -0.5 * ((Beta + Gamma) - sqrt((Beta - Gamma).^2 + 4*C));
% 
% % Plot the surface
% figure;
% surf(Beta, Gamma, Rate);
% xlabel('\beta');
% ylabel('\gamma');
% zlabel('Rate');
% title('Theoretical Rate Surface with C = 1');
% shading interp;      % smoother surface
% colorbar;
% 
% %%
% gamma_vals = linspace(0, 5, 200);
% 
% % Choose some representative beta values
% beta_list = [0.5, 1, 2, 3, 4, 5];
% 
% % Fixed constant
% C = 1;
% 
% % Initialize figure
% figure;
% hold on;
% 
% % Plot rate vs gamma for each beta
% for i = 1:length(beta_list)
%     beta = beta_list(i);
%     rate_vals = get_theoretical_rate_exponential(beta, gamma_vals, C);
%     plot(gamma_vals, rate_vals, 'LineWidth', 1.5, ...
%          'DisplayName', sprintf('\\beta = %.1f', beta));
% end
% 
% xlabel('\gamma');
% ylabel('Rate');
% title('Rate vs. \gamma for Various \beta Values (C = 1)');
% legend('Location', 'best');
% grid on;