function [fval, grad] = usrfntry(par, C, lambda1, beta1, beta2)
% par = data.par;
% C = data.C;
% lambda1 = data.lambda1;
% beta1 = data.beta1;
% beta2 = data.beta2;
p = C' * par;
M = abs(p) - lambda1./beta1;
M(M < 0) = 0;
fval = beta2' * M.^2;
grad = 2*C * (M.*sign(p).*beta2);
fval = [fval;0;0];
end