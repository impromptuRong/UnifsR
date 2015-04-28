function [model, score, alpha] = unifsR_con(TRAIN, y, TEST, beta1, beta2, lambda1, D, x0)
%% Solve minimum convex-hull distance
[npar, nfea] = size(TRAIN);
noZero = (full(sum(TRAIN,1))~=0) & (beta1'~=0);
C = diag(y)*TRAIN(:,noZero);
beta1_ori = beta1;
beta1 = beta1(noZero);
beta2 = beta2(noZero);
lb = zeros(npar,1);
ub = [y<0,y>0]*D*100;
A_eq = [ones(1,npar); y'];
b_eq = [200; 0];
alpha = 1;

%% Optimization Formula %%
% sum(bw*(rho-(alpha*R)./bw).^2) = lambda^2*(alpha*R)^2
% min sum([rho(k)-lambda1/beta1(k)](+)^2*beta2(k)), lambda1=alpha*R, fval=(R-alpha*R)^2
    function [fval, grad] = fn(par)
        p = C' * par;
        M = abs(p) - lambda1./beta1;
        M(M < 0) = 0;
        fval = beta2' * M.^2;
        grad = 2*C * (M.*sign(p).*beta2);
    end
    function fval = zf(par)
        [fval, ~] = fn(par);
    end
    function grad = zg(par)
        [~, grad] = fn(par);
    end
%     function [fval, grad] = fnn(par)
%         [fval, grad] = fn(par);
%         fval = [fval;0;0];
%     end

%% Solve quadratic system for lambda1=0 (alpha=0) and generating initial point
if(~isempty(lambda1))
    if(lambda1 == 0 || isempty(x0))
        H = C*diag(beta2)*C';
%         options = optimset('Algorithm','interior-point-convex','Display','off');
%         [u, R, ~] = quadprog(H+H', [], [], [], A_eq, b_eq, lb, ub, [], options);
        Prob = qpAssign(H+H', zeros(npar,1), A_eq, b_eq, b_eq, lb, ub, [], 'QP');
        Result = tomRun('qpopt', Prob, 0);
        u = Result.x_k;
        R = Result.f_k;
        x0 = u;
    end;
    
    if(lambda1 > 0)
%         options = optimset('Algorithm','sqp','Display','off','GradObj','on','UseParallel','always');
%         [u, R, ~] = fmincon(@fn, x0, [], [], A_eq, b_eq, lb, ub, [], options);
%         [u1, R1, ~] = snopt(x0,lb,ub,[0;b_eq],[Inf;b_eq],'fnn',[ones(npar,1);y],[ones(npar,1)*2;ones(npar,1)*3],[1:npar,1:npar]',ones(npar,1),(1:npar)');
        Prob = conAssign(@zf, @zg, [], [], lb, ub, 'NLP', x0, [], [], A_eq, b_eq, b_eq);
        Result = tomRun('npsol', Prob, 0);
        u = Result.x_k;
        R = Result.f_k;
    end
    R = lambda1 + sqrt(R/nfea);
    alpha = lambda1/R;
end

%% Solve linear system for lambda1=max (alpha=1) and generating initial point if not provided
if(alpha > 0.999)
    [model, score, alpha] = unifsR_lin(TRAIN, y, TEST, beta1_ori, D);
    return;
end

%% Calculate Support Vector, separation hyperplane
tol = 1e-8;
SV = zeros(1,npar);
SV(u < tol) = -1;
SV(u > ub-tol) = 1;

p = C' * u;
rw = abs(p) - R*alpha./beta1;
rw(rw < 0) = 0;
rw = rw.*sign(p).*beta2;  % sum(w.^2./beta2/N) - lambda2^2;
rw = rw/(rw'*p)*R;

A_mat = [(y+1)/2, (y-1)/2];
b_mat = full(C * rw);
A_ineq = [A_mat(SV==-1,:); -A_mat(SV~=-1,:)];
b_ineq = [b_mat(SV==-1); -b_mat(SV~=-1)];
% min(b_ineq(A_ineq(:,1)==1)) >= a >= max(-b_ineq(A_ineq(:,1)==-1)); min(b_ineq(A_ineq(:,2)==-1)) >= -b >= max(-b_ineq(A_ineq(:,2)==1));
opt_var = [max(-b_ineq(A_ineq(:,1)==-1)); min(b_ineq(A_ineq(:,2)==1))];
dsplane = 100 * (opt_var(1)-opt_var(2));
rw = [-opt_var(1)-opt_var(2);2*rw];
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
