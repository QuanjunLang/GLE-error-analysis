% function [f1, f2, K1, K2] = switch_kernel(type)
%     switch type
%         case 'exp'
%             tau1 = 3;
%             tau2 = 2.8;
%             f1 = @(t) -exp(-t * tau1);
%             f2 = @(t) -exp(-t * tau2);
%             % K1 = @(t, s) -exp(-(t - s) * tau1) .* (t >= s);
%             % K2 = @(t, s) -exp(-(t - s) * tau2) .* (t >= s);
%         otherwise
%             error('Unknown kernel type: %s', type);
%     end
%     K1 = @(t, s) f1(t-s).* (t >= s);
%     K2 = @(t, s) f2(t-s).* (t >= s);
% end
% 

function [f1, f2, K1, K2] = switch_kernel(type)
    switch type
        case 'exp'
            tau1 = 1; tau2 = 1.2;
            f1 = @(tau) -exp(-tau * tau1);
            f2 = @(tau) -exp(-tau * tau2);
        case 'power'
            alpha1 = 3; alpha2 = 3.1;
            f1 = @(tau) -(1 + tau).^(-alpha1);
            f2 = @(tau) -(1 + tau).^(-alpha2);
        case 'osc'
            f1 = @(tau) exp(-10*tau) .* sin(tau).* (tau <= 10);
            f2 = @(tau) 10 * exp(-10*tau) .* sin(tau).* (tau <= 10);
        otherwise
            error('Unknown kernel type');
    end
    K1 = @(t, s) f1(t - s) .* (t >= s);
    K2 = @(t, s) f2(t - s) .* (t >= s);
end