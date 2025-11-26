function result = ensemble_test_S_c2_C2(sysInfo, f2_info)


kernel_type = 'powerlaw_test';                       % Specify the type of kernel to generate ('osc', 'exp', 'powerlaw', 'diff')
decay = 4;
amp1 = 1;
amp2 = 1;
intercept1 = 1;
intercept2 = 1;



kernelInfo = new_kernel_info(kernel_type, sysInfo, f2_info, 'decay', decay, 'amp1', amp1, 'amp2', amp2, 'intercept1', intercept1, 'intercept2', intercept2); % Generate the kernel functions
% plot_kernels(kernelInfo);                       % Plot the kernel functions in three panels: linear, semilog-y, and log-log


threshold = get_theoretical_threshold_powerlaw(decay, amp1, intercept1);
fprintf('For powerlaw kernel, the threshold for gamma is : %.2f\n', threshold)


%%


f_diff = @(t) abs(kernelInfo.f1(t) - kernelInfo.f2(t));
h = @(t) (t + intercept1).^(-2*decay+2);
S_norm_delta_K = Schur_norm(f_diff, h);
S_norm_tilde_K = Schur_norm(kernelInfo.f2, h);


%% Generate Ensemble trajectories
M = 20;
traj_1 = zeros(M, sysInfo.d, sysInfo.tn);
traj_2 = zeros(M, sysInfo.d, sysInfo.tn);

v0_sigma = 0.5;
v0_mu = 1;

for i = 1:M

    noise = get_noise(sysInfo);

    traj_input{1}.v0 = v0_mu + randn(sysInfo.d, 1)*v0_sigma;
    traj_input{1}.noise = noise;

    temp = generate_data_first_order_highdim(sysInfo, kernelInfo.K1, traj_input{1});
    traj_1(i, :, :) = temp.v;

    % traj_input{2}.v0 = v0_mu + randn(sysInfo.d, 1)*v0_sigma;
    traj_input{2}.v0 = traj_input{1}.v0;
    traj_input{2}.noise = noise;
    temp = generate_data_first_order_highdim(sysInfo, kernelInfo.K2, traj_input{2});
    traj_2(i, :, :) = temp.v;
end

traj_1_mean = squeeze(mean(traj_1.^2, 1));
traj_2_mean = squeeze(mean(traj_2.^2, 1));
traj_diff = squeeze(mean((traj_1 - traj_2).^2, 1));

%%
% fitOptions.fit_t0 = 4;
% fitOptions.fit_v0 = 1e-10;
% fitOptions.plotON = true;
% fitOptions.type = 'loglog';    % semilog (Exp) or loglog (Powerlaw) or all
% fitOptions.displayOn = true;


% fitOptions_1        = fitOptions;
% fitOptions_1.ttl    = 'traj 1';
% fitOptions_1.plotON = false;
% fitInfo_1 = fit_traj_rate(traj_1_mean, sysInfo, fitOptions_1);
% 
% fitOptions_2        = fitOptions;
% fitOptions_2.ttl    = 'traj 2';
% fitOptions_2.plotON = false;
% fitInfo_2 = fit_traj_rate(traj_2_mean, sysInfo, fitOptions_2);
% 
% 

%%
% fitOptions_diff        = fitOptions;
% fitOptions_diff.ttl    = 'traj difference';
% fitOptions_diff.fit_v0 = 1e-8;
% fitOptions_diff.fit_t0 = 2;
% fitOptions_diff.h = h;
% fitOptions_diff.displayOn = false;
% fitOptions_diff.plotON = false;
% 
% 
% fitInfo_diff = fit_traj_rate(traj_diff, sysInfo, fitOptions_diff);


%%
d = sysInfo.d;
if d  == 1
    vnorm = abs(traj_diff)';
else

    vnorm = vecnorm(traj_diff, 2, 1);  % L2 norm at each time
end


ht = h(sysInfo.t_grid);

range_c2 = sysInfo.tn - 50;
% c2 = mean(vnorm(range_c2:end));
c2 = sysInfo.sigma.^2;
% 
% c2 = 0;
C2 = max(vnorm./(ht + c2));
%%

fig = figure;


subplot(121);
plot(sysInfo.t_grid, vnorm./(ht + c2))







subplot(122);hold on;grid on;
loglog(sysInfo.t_grid, ht, 'DisplayName','h');

loglog(sysInfo.t_grid, C2*(ht + c2), 'DisplayName','upper bound')
loglog(sysInfo.t_grid, vnorm, 'DisplayName','E|W_t^2');



% C_loglog = exp(fitInfo_diff.logC_loglog);
% alpha_loglog = fitInfo_diff.alpha_loglog;
% vnorm_fitline_loglog = C_loglog * t.^alpha_loglog;      % Fitted line: vnorm ≈ C * t^alpha
% loglog(t, vnorm_fitline_loglog, 'm-', 'LineWidth', 1.5, 'DisplayName','Fitted rate');

set(gca, 'YScale', 'log')
legend()

% close(fig);

%%
result.C2 = C2;
result.c2 = c2;
result.S_norm_delta_K = S_norm_delta_K;
result.S_norm_tilde_K = S_norm_tilde_K;
result.vnorm = vnorm;

end




