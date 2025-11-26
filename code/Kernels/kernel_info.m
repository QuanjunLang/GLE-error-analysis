function kernelInfo = kernel_info(kernel_type, sysInfo, varargin)
    % Create parser
    p = inputParser;
    validKernelTypes = {'osc', 'exp', 'powerlaw', 'diff', 'powerlaw_sharp', 'powerlaw_test'};
    
    % Required positional input
    addRequired(p, 'kernel_type', @(x) any(validatestring(x, validKernelTypes)));
    addRequired(p, 'sysInfo');
    % Optional name-value parameters
    addParameter(p, 'decay', 2, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'amp1', 1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'amp2', 1.1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'intercept1', 1, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'intercept2', 1, @(x) isnumeric(x) && isscalar(x));
    
    
    % Parse inputs
    parse(p, kernel_type, sysInfo, varargin{:});
    kernel_type = p.Results.kernel_type;
    sysInfo = p.Results.sysInfo;
    decay = p.Results.decay;
    amp1 = p.Results.amp1;
    amp2 = p.Results.amp2;

    intercept1 = p.Results.intercept1;
    intercept2 = p.Results.intercept2;

    % Define kernels based on type
    switch kernel_type
        case "osc"
            f1 = @(tau) amp1 * exp(-decay * tau) .* sin(0.5 * tau + pi);
            % f2 = @(tau) -amp2 * exp(-decay*2 * tau) .* sin(0.3 * tau + 2) .*(tau < 3);
            f2 = @(tau) amp2 * exp(-decay*2 * tau) .* sin(0.3 * tau + 0.3);
            label1 = sprintf('$f_1(t) = -%.2f e^{-%dt} \\sin(0.5t+1)$', amp1, decay);
            label2 = sprintf('$f_2(t) = -%.2f e^{-%dt} \\sin(0.5t+1)$', amp2, decay);
            title_str = "Oscillatory Exponential Kernels";
        case "exp"
            f1 = @(tau) -amp1 * exp(-decay * tau);
            f2 = @(tau) -amp2 * exp(-decay * tau);
            label1 = sprintf('$f_1(t) = -%.2f e^{-%dt}$', amp1, decay);
            label2 = sprintf('$f_2(t) = -%.2f e^{-%dt}$', amp2, decay);
            title_str = "Exponential Kernels";

        case "powerlaw"
            f1 = @(tau) -amp1 * (tau + intercept1).^(-decay);
            f2 = @(tau) -amp2 * (tau + intercept2).^(-decay);
            label1 = sprintf('$f_1(t) = -%.2f(t+%.2f)^{-%d}$', amp1, intercept1, decay);
            label2 = sprintf('$f_2(t) = -%.2f(t+%.2f)^{-%d}$', amp2, intercept2, decay);
            title_str = "Power-law Kernels";

        case "diff"
            f1 = @(tau) amp1 * (tau + 3).^(-0.5);
            f2 = @(tau) amp1 * exp(-0.1 * tau) .* sin(tau);
            label1 = sprintf('$f_1(t) = %.2f(t+3)^{-0.5}$', amp1);
            label2 = sprintf('$f_2(t) = %.2f e^{-0.1t} \\sin(t)$', amp1);
            title_str = "Different Type Kernels";

        case "powerlaw_osc"
            f1 = @(tau) -amp1 * (tau + intercept1).^(-decay).* sin(0.5 * tau + pi);
            f2 = @(tau) -amp2 * (tau + intercept2).^(-decay-0.5).* sin(0.5 * tau + pi);
            label1 = sprintf('$f_1(t) = -%.2f(t+%.2f)^{-%d}$', amp1, intercept1, decay);
            label2 = sprintf('$f_2(t) = -%.2f(t+%.2f)^{-%d}$', amp2, intercept2, decay);
            title_str = "Power-law Kernels";
        case "powerlaw_test"
            f1 = @(tau) amp1 * (tau + intercept1).^(-decay);
            f2 = @(tau) amp2 * (tau + intercept2).^(-decay).*(tau < 2);
            label1 = sprintf('$f_1(t) = -%.2f(t+%.2f)^{-%d}$', amp1, intercept1, decay);
            label2 = sprintf('$f_2(t) = -%.2f(t+%.2f)^{-%d}$', amp2, intercept2, decay);
            title_str = "Power-law Kernels";
    end

    d = sysInfo.d;
    rng(1);
    C = randn(d, d);
    A = C*C';
    [U, ~, ~] = eig(A);

    % Define full kernels K(t, s)
    K1 = @(t, s) U*f1(t-s)*U';
    K2 = @(t, s) U*f2(t-s)*U';

    kernelInfo.f1 = f1;
    kernelInfo.f2 = f2;
    kernelInfo.K1 = K1;
    kernelInfo.K2 = K2;


    kernelInfo.K1_norm = @(t, s) norm(K1(t, s));
    kernelInfo.K2_norm = @(t, s) norm(K2(t, s));

    
    kernelInfo.label1 = label1;
    kernelInfo.label2 = label2;
    kernelInfo.title_str = title_str;

    % kernelInfo.f_diff = @(t) abs(f1(t) - f2(t));
end