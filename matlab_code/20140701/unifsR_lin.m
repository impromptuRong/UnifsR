function [model, score, alpha] = unifsR_lin(TRAIN, y, TEST, beta1, D)
%% C, 100*u, 100*R, 100*d, 1*w
[npar, nfea] = size(TRAIN);
noZero = (full(sum(TRAIN,1))~=0) & (beta1'~=0);
C = diag(y)*TRAIN(:,noZero);
beta1 = beta1(noZero);
np = sum(noZero);

%% Solve linear system for lambda1=max (alpha=1) and generating initial point if not provided
lb = zeros(npar+1,1);
ub = [y<0,y>0]*D*100;
A_eq = [ones(1,npar),0; y',0];
b_eq = [200; 0];
A_ineq = [[diag(beta1)*C'; -diag(beta1)*C'],-ones(2*np,1)];
b_ineq = zeros(2*np,1);

%% Matlab linprog linear solver %%
% options = optimset('Display','off','Algorithm','Active-set');
% [u, ~, ~] = linprog([zeros(npar,1);1], A_ineq, b_ineq, A_eq, b_eq, lb, [ub;Inf], [], options);
%% TOMLAB linear programming %%
Prob = lpAssign(sparse(npar+1,1,1), [A_ineq; A_eq], [-Inf(size(b_ineq));b_eq], [b_ineq;b_eq], lb, [ub;Inf], [], 'linear');
Result = tomRun('lpopt', Prob, 0);
u = Result.x_k;
%% snopt solver %%
% uf = @(x) A*x;
% A = [zeros(1,npar),1;A_ineq;A_eq];
% [jAvar, iAfun] = meshgrid(1:size(A,2),1:size(A,1));
% [u, ~, ~] = snopt(zeros(npar+1,1), lb, [ub;Inf], [0;-Inf(size(b_ineq));b_eq], [Inf;b_ineq;b_eq], 'uf', A(:), iAfun(:), jAvar(:), [], []);

%% Calculate Support Vector, separation hyperplane
alpha = 1;
R = u(end);
u = u(1:npar);
tol = 1e-6;
SV = zeros(1,npar);
SV(u < tol) = -1;
SV(u > ub-tol) = 1;
p = C' * u;
rescale = abs(p) > R./beta1 - tol*alpha;
nres = sum(rescale);
ps = sign(p(rescale));

A_mat = full([(y+1)/2, (y-1)/2, -C(:,rescale)]);
b_mat = zeros(npar,1);
reduce = rref([A_mat(SV==0,:),b_mat(SV==0)], 1e-2);
reduce = reduce(sum(reduce(:,1:end-1),2)~=0,:);
A_eq = [0,0,ps'./beta1(rescale)';reduce(:,1:end-1)];
b_eq = [1;reduce(:,end)];
A_ineq = [A_mat(SV==-1,:); -A_mat(SV==1,:)];
b_ineq = [b_mat(SV==-1); -b_mat(SV==1)];
tmp = [[-Inf;-Inf;zeros(nres,1)], [Inf;Inf;ps]];
lb = min(tmp,[],2);
ub = max(tmp,[],2);

%% Matlab linprog linear solver %%
% options = optimset('Display','off','Algorithm','Active-set');
% [opt_var, dsplane, ~] = linprog([100;-100;zeros(nres,1)],A_ineq,b_ineq,A_eq,b_eq,lb,ub,[],options);
%% TOMLAB linear programming %%
Prob = lpAssign([100;-100;zeros(nres,1)], [A_ineq;A_eq], [-Inf(size(b_ineq));b_eq], [b_ineq;b_eq], lb, ub, [], 'linear');
Result = tomRun('lpopt', Prob, 0);
opt_var = Result.x_k;
dsplane = Result.f_k;
%% snopt solver %%
% A = [[100,-100,zeros(1,nres)];A_ineq;A_eq];
% [jAvar, iAfun] = meshgrid(1:size(A,2),1:size(A,1));
% [opt_var, dsplane, ~] = snopt(zeros(nres+2,1), lb, ub, [0;-Inf(size(b_ineq));b_eq], [Inf;b_ineq;b_eq], 'uf', A(:), iAfun(:), jAvar(:), [], []);
% dsplane = dsplane(1);

%% Calculating w and generate model %%
rw = [-opt_var(1)-opt_var(2);zeros(np,1)];
rw([false;rescale]) = 2*opt_var(3:end);
w = zeros(nfea+1,1);
w([true,noZero]) = rw/dsplane*100;

model.u = u;
model.R = R;
model.w = w;
model.d = dsplane;
model.D = D;
model.alpha = alpha;
score.train = [ones(npar,1),TRAIN]*w + 10^-20;
score.test = [ones(size(TEST,1),1),TEST]*w + 10^-20;
end
