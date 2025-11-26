function result = ensemble_test_second_order(sysInfo, perturbInfo, kernelInfo, potentialInfo, LyapunovInfo, fitInfo)


%% update the kernel info and compute kernel error
kernelInfo = new_kernel_info_structure(kernelInfo, sysInfo, perturbInfo);


%% Generate Ensemble trajectories
M = sysInfo.M;
trajInfo = cell(2, M);
v0_sigma = 0.5;
v0_mu = 1;

r_square_diff   = zeros(sysInfo.tn, M);
x_diff          = zeros(sysInfo.d, sysInfo.tn, M);
v_diff          = zeros(sysInfo.d, sysInfo.tn, M);


for i = 1:M
    % Generate trajectory
    fprintf('Generating the trajectory: No.%d\n', i)
    noise = get_noise(sysInfo);

    traj_input{1}.x0 = v0_mu + randn(sysInfo.d, 1)*v0_sigma;
    traj_input{1}.v0 = v0_mu + randn(sysInfo.d, 1)*v0_sigma;
    traj_input{1}.noise = noise;
    trajInfo{1, i} = generate_data_second_order_highdim(sysInfo, kernelInfo.K1, potentialInfo, traj_input{1});

    traj_input{2}.x0 = traj_input{1}.x0;
    traj_input{2}.v0 = traj_input{1}.v0;
    traj_input{2}.noise = noise;
    trajInfo{2, i} = generate_data_second_order_highdim(sysInfo, kernelInfo.K2, potentialInfo, traj_input{2});
    
    r_square_diff(:, i) = compute_Lyapunov_distance(trajInfo{1, i}.x, trajInfo{1, i}.v, trajInfo{2, i}.x, trajInfo{2, i}.v, LyapunovInfo);
    x_diff(:, :, i) = trajInfo{1, i}.x - trajInfo{2, i}.x;
    v_diff(:, :, i) = trajInfo{1, i}.v - trajInfo{2, i}.v;
end


%% Fit and find c2 C2
xt = sqrt(mean(x_diff.^2, [1, 3]))';
vt = sqrt(mean(x_diff.^2, [1, 3]))';

Rt = mean(r_square_diff, 2);
ht = kernelInfo.h(sysInfo.t_grid)';

c2 = sysInfo.sigma.^2;
C2 = max(Rt./(ht + c2));
%%


% fitInfo.t0 = 10;
% fitInfo.v0 = sysInfo.sigma^3+eps;
% fitInfo.plotON = false;
% fitInfo.printON = true;

fitInfo.Rt      = Rt;
fitInfo.ht      = ht;
fitInfo.c2      = c2;
fitInfo.C2      = C2;



% fitON = ~(C2 < 1e-5);
fitON = ~(max(Rt)<fitInfo.v0);

if fitON
    fit_result = fit_new_exp(sysInfo, fitInfo);
    result.fit_result = fit_result;
else
    result.fit_result = [];
end


%%
result.C2 = C2;
result.c2 = c2;
result.S_norm_delta_K = kernelInfo.S_norm_delta_K;
result.S_norm_tilde_K = kernelInfo.S_norm_tilde_K;
result.Rt = Rt;




%%



end




