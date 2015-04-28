function [model, score, convexD] = unifsR_init0(TRAIN, y, TEST, beta2, C)
%% C, 100*u, 100*R, 100*d, 1*w
[npar, nfea] = size(TRAIN);
noZero = (full(sum(TRAIN,1))~=0) & (beta2'~=0);
E = diag(y)*TRAIN(:,noZero);
beta2 = beta2(noZero);

%% Solve quadratic system for lambda1=0 (alpha=0) and generating initial point
H = E*diag(beta2)*E';
lC = zeros(npar,1);
uC = C*ones(npar,1);
% options = optimset('Algorithm','interior-point-convex','Display','off');
% [u, ~, ~] = quadprog(H+H', -2*ones(1,npar), [], [], y', 0, lC, uC, [], options);
Prob = qpAssign(H+H', -2*ones(npar,1), y', 0, 0, lC, uC, [], 'QP');
Result = tomRun('qpopt', Prob, 0);
u = Result.x_k;
rw = E'*u.*beta2;

%% Calculate Support Vector, separation hyperplane
tol = 1e-8;
SV = zeros(1,npar);
SV(u < tol) = -1;
SV(u > uC-tol) = 1;

if(any(~SV))
    b = mean((1-E(SV==0,:)*rw).*y(SV==0));
else
    A_ineq = sign((SV'+0.1).*y);
    b_ineq = sign(SV'+0.1).*(1 - E * rw);
    b = (min(b_ineq(A_ineq==1)) + max(-b_ineq(A_ineq==-1)))/2;
end
w = zeros(nfea+1,1);
w([true,noZero]) = [b;rw];

delta = sum(u)/200;
convexD = uC/delta;
u = u/delta;
R = sqrt(u'*H*u/nfea);
dsplane = 200/sqrt(nfea*sum(rw.*rw./beta2));
alpha = 0;
D = max(u)/100;

%% Output model and scores %%
model.u = u;
model.R = R;
model.w = w;
model.d = dsplane;
model.D = [D,max(convexD)];
model.delta = delta;
model.alpha = alpha;
score.train = [ones(npar,1),TRAIN]*w + 10^-20;
score.test = score.train;
if(~isempty(TEST))
    score.test = [ones(size(TEST,1),1),TEST]*w + 10^-20;
end
end