function kernelInfo = new_kernel_info_structure(kernelInfo, sysInfo, f2_info)


d           = sysInfo.d;
type        = kernelInfo.type;
amp         = kernelInfo.amp;
decay       = kernelInfo.decay;
intercept   = kernelInfo.intercept;

%%
s = rng;
rng(42); 
C = randn(d, d);
A = C*C';
[U, S, ~] = eig(A);
rng(s);

S_diag = diag(S) + (decay - min(diag(S)));         % This is to make sure that the lowest spectrum is the decay parameter

S_diag = decay*ones(d, 1);
% Here we only consider the case of convolution kernel
% k is the convolution kernel, and K(t, s) = k(t-s)
% f is the 1D example of showing the perturbation effect

switch type
    case "exp_test"
        k1 = @(tau) amp * U * diag(exp(-S_diag * tau)) * U';
        switch f2_info.method
            case 'shift'
                k2 = @(tau) amp* U * diag(exp(-(S_diag + f2_info.para) * tau)) * U';
            case 'rate_plus'
                k2 = @(tau) amp* U * diag(exp(-S_diag * (tau + f2_info.para))) * U';
            case 'cut'
                k2 = @(tau) amp* U * diag(exp(-S_diag * tau)) .*(tau < f2_info.para) * U';
            case 'osc'
                k2 = @(tau) amp* U * diag(exp(-S_diag * tau)) .* cos(tau * f2_info.para) * U';
        end

        f1 = @(tau) amp * exp(- decay * tau);
        switch f2_info.method                   % Example showing the effect of the perturbation.
            case 'shift'
                f2 = @(tau) amp* exp(-(decay + f2_info.para) * tau);
            case 'rate_plus'
                f2 = @(tau) amp* exp(-decay * (tau + f2_info.para));
            case 'cut'
                f2 = @(tau) amp* exp(-decay * tau) .*(tau < f2_info.para);
            case 'osc'
                f2 = @(tau) amp* exp(-decay * tau) .* cos(tau * f2_info.para);
        end


        label1 = '$f_1(t)$';
        label2 = '$f_2(t)$';
        title_str = "Exponential Kernels";
        K1 = @(t, s) k1(t-s);
        K2 = @(t, s) k2(t-s);

    case "powerlaw_test"
        k1 = @(tau) amp * U * diag((tau + intercept).^(-S_diag)) * U';
        switch f2_info.method
            case 'shift'
                k2 = @(tau) amp * U * diag((tau + intercept + f2_info.para).^(-S_diag)) * U';
            case 'rate_plus'
                k2 = @(tau) amp * U * diag((tau + intercept).^(-S_diag - f2_info.para)) * U';
            case 'cut'
                k2 = @(tau) amp * U * diag((tau + intercept).^(-S_diag).*(tau < f2_info.para)) * U';
            case 'osc'
                k2 = @(tau) amp * U * diag((tau + intercept).^(-S_diag).* cos(tau * f2_info.para)) * U';
        end

        f1 = @(tau) amp * (tau + intercept).^(-decay);
        switch f2_info.method
            case 'shift'
                f2 = @(tau) amp * (tau + intercept + f2_info.para).^(-decay);
            case 'rate_plus'
                f2 = @(tau) amp * (tau + intercept).^(-decay - f2_info.para);
            case 'cut'
                f2 = @(tau) amp * (tau + intercept).^(-decay).*(tau < f2_info.para);
            case 'osc'
                f2 = @(tau) amp * (tau + intercept).^(-decay) .* cos(tau * f2_info.para);
        end

        label1 = '$f_1(t)$';
        label2 = '$f_2(t)$';
        title_str = "Power-law Kernels";
        K1 = @(t, s) k1(t-s);
        K2 = @(t, s) k2(t-s);
end


%%

kernelInfo.f1 = f1;
kernelInfo.f2 = f2;
kernelInfo.K1 = K1;
kernelInfo.K2 = K2;


kernelInfo.K1_norm = @(t, s) norm(K1(t, s));
kernelInfo.K2_norm = @(t, s) norm(K2(t, s));


kernelInfo.label1 = label1;
kernelInfo.label2 = label2;
kernelInfo.title_str = title_str;

%% Compute kernel erorr using smalle k1 k2
% In the case of matrix exponential, the matrix operator norm is precicesly
% just the norm using f1 and f2
h = kernelInfo.h;
kernelInfo.S_norm_true_K = Schur_norm(f1, h);
kernelInfo.S_norm_tilde_K = Schur_norm(f2, h);
kernelInfo.S_norm_delta_K = Schur_norm(@(t) f1(t) - f2(t), h);
end