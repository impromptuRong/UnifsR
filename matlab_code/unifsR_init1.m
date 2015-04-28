function [model, score, convexD] = unifsR_init1(TRAIN, y, TEST, beta1, C)
%% C, 100*u, 100*R, 100*d, 1*w
[npar, nfea] = size(TRAIN);
noZero = (full(sum(TRAIN,1))~=0) & (beta1'~=0);
E = diag(y)*TRAIN(:,noZero);
beta1 = beta1(noZero);
np = sum(noZero);

%% Solve linear system for lambda1=max (alpha=1) and generating initial point if not provided
lC = zeros(npar,1);
uC = C*ones(npar,1);
A_eq = y';
b_eq = 0;
A_ineq = diag(beta1)*E';
b_ineq = ones(np,1);

%% Matlab linprog linear solver %%
% options = optimset('Display','off','Algorithm','Active-set');
% [u, ~, ~] = linprog(-ones(npar,1), [A_ineq;-A_ineq], [b_ineq;b_ineq], A_eq, b_eq, lC, uC, [], options);
%% TOMLAB linear programming %%
Prob = lpAssign(-ones(npar,1), [A_ineq;A_eq], [-b_ineq;b_eq], [b_ineq;b_eq], lC, uC, [], 'linear');
Result = tomRun('lpopt', Prob, 0);
u = Result.x_k;
%% snopt solver %%
% uf = @(x) A*x;
% A = [zeros(1,npar),1;A_ineq;A_eq];
% [jAvar, iAfun] = meshgrid(1:size(A,2),1:size(A,1));
% [u, ~, ~] = snopt(zeros(npar+1,1), lC, [uC;Inf], [0;-Inf(size(b_ineq));b_eq], [Inf;b_ineq;b_eq], 'uf', A(:), iAfun(:), jAvar(:), [], []);

%% Calculate Support Vector, separation hyperplane
tol = 1e-8;
SV = zeros(1,npar);
SV(u < tol) = -1;
SV(u > uC-tol) = 1;
p = A_ineq*u;
rescale = abs(p) > 1-tol;
nres = sum(rescale);

A_mat = [y,E(:,rescale)];
b_mat = ones(npar,1);
reduce = rref([A_mat(SV==0,:), b_mat(SV==0)], 1e-2);
reduce = reduce(sum(reduce(:,1:end-1),2)~=0,:);
A_eq = reduce(:,1:end-1);
b_eq = reduce(:,end);
A_ineq = [-A_mat(SV==-1,:); A_mat(SV==1,:)];
b_ineq = [-b_mat(SV==-1); b_mat(SV==1)];
tmp = [-Inf,Inf;zeros(nres,1),Inf*sign(p(rescale))];
lb = min(tmp,[],2);
ub = max(tmp,[],2);

%% Matlab linprog linear solver %%
% options = optimset('Display','off','Algorithm','Active-set');
% [opt_var, ~, ~] = linprog([0;sign(p(rescale))],A_ineq,b_ineq,A_eq,b_eq,lb,ub,[],options);
%% TOMLAB linear programming %%
Prob = lpAssign([0;sign(p(rescale))], [A_ineq;A_eq], [-Inf(size(b_ineq));b_eq], [b_ineq;b_eq], lb, ub, [], 'linear');
Result = tomRun('lpopt', Prob, 0);
opt_var = Result.x_k;
%% snopt solver %%
% A = [[100,-100,zeros(1,nres)];A_ineq;A_eq];
% [jAvar, iAfun] = meshgrid(1:size(A,2),1:size(A,1));
% [opt_var, dsplane, ~] = snopt(zeros(nres+2,1), lb, ub, [0;-Inf(size(b_ineq));b_eq], [Inf;b_ineq;b_eq], 'uf', A(:), iAfun(:), jAvar(:), [], []);
% dsplane = dsplane(1);

%% Calculating w and generate model %%
rw = zeros(np+1,1);
rw([true;rescale]) = opt_var;
w = zeros(nfea+1,1);
w([true,noZero]) = rw;

delta = sum(u)/200;
convexD = uC/delta;
u = u/delta;
R = max(abs(diag(beta1)*E'*u));
dsplane = 200/sum(abs(opt_var(2:end))./beta1(rescale));
alpha = 1;
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