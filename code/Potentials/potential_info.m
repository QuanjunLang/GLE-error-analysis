function potentialInfo = potential_info(sysInfo, potentialInfo)

x_min_temp  = potentialInfo.x_min_temp;
u           = potentialInfo.u;
kappa       = potentialInfo.kappa;
type        = potentialInfo.type;
d           = sysInfo.d;

switch type
    case 'smoothed_huber'
        delta = 5;

        C = rand(d, d) + 1;
        [V, ~, ~] = eig(C);
        R = V*(diag(rand(d, 1)+1))*V';
        kappa = min(eig(R));


        G       = @(x) sum(sqrt(x.^2 + delta^2) - delta); % Potential function
        G_grad  = @(x) x ./ sqrt(x.^2 + delta^2);    % Gradient of G

        U       = @(x) x' * R * x/2 + G(x-x_min_temp);
        U_grad  = @(x) R*x + G_grad(x-x_min_temp);

    case 'quadratic'
        delta = 5;

        R = eye(d);
        kappa = 1;

        G       = @(x) sum(sqrt(x.^2 + delta^2) - delta); % Potential function
        G_grad  = @(x) x ./ sqrt(x.^2 + delta^2);    % Gradient of G

        U       = @(x) x' * R * x/2;
        U_grad  = @(x) R*x;

    case 'quadratic_new'
        delta = 100;

        A = rand(d, d);
        R = (A*A');
        [UU, S] = eig(R);
        S = diag(diag(S) + (kappa - min(diag(S)))); 
        R = UU*S*UU';

        G       = @(x) sum(sqrt((x-x_min_temp).^2 + delta^2) - delta); % Potential function
        G_grad  = @(x) (x-x_min_temp) ./ sqrt((x-x_min_temp).^2 + delta^2);    % Gradient of G

        U       = @(x) x' * R * x/2 + G(x);
        U_grad  = @(x) R*x + G_grad(x);

        LG = 1/delta;
end


% Optimization
fun = @(x) deal(U(x), U_grad(x));
x0 = zeros(size(x_min_temp));
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'GradObj', 'on', 'Display', 'off');
[x_opt, ~] = fminunc(fun, x0, options);


potentialInfo.u     = u;
potentialInfo.LG    = LG;
potentialInfo.kappa = kappa;
potentialInfo.R     = R;

potentialInfo.U         = U;
potentialInfo.U_grad    = U_grad;
potentialInfo.G         = G;
potentialInfo.G_grad    = G_grad;


potentialInfo.U_argmin = x_opt;
end


