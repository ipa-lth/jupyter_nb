function [k0, k1] = k_explizit_Ab3(roots_p1, roots_pmin, pmin, A, b)

n = size(roots_p1, 2);

# named 'a_start' at Boris:calcCntrlRLStyle.m
a_tilde_p1 = fliplr(poly(roots_p1)(2:n+1)); # get the characteristical polynomial backwards

# named 'a_stop' at Boris:calcCntrlRLStyle.m
a_tilde_pmin = fliplr(poly(roots_pmin)(2:n+1)); # get the characteristical polynomial backwards

M = [ kron(eye(n),[1 1]') , kron(eye(n),[1 1/pmin]') ];
a_p1_pmin = kron(a_tilde_p1',[1 0]') + kron(a_tilde_pmin',[0 1]');

a_hat = M \ a_p1_pmin;
a_hat_2 = a_hat(1:end/2);
a_hat_1 = a_hat(end/2+1:end);

[A_R, _, _, _ , T, _] = get_Steuerungsnormalform(A, b, b, 0);
a = A_R(size(A_R, 1), :)';

k1 = a_hat_1;
k0 = a_hat_2 + a;

k1 = T'*k1;
k0 = T'*k0;

return