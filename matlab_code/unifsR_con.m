function [model, score, alpha] = unifsR_con(TRAIN, y, TEST, beta1, beta2, lambda1, convexD, x0)
%% Solve minimum convex-hull distance
[npar, nfea] = size(TRAIN);
noZero = (full(sum(TRAIN,1))~=0) & (beta1'~=0) & (beta2'~=0);
E = diag(y)*TRAIN(:,noZero);
beta1_ori = beta1;
beta1 = beta1(noZero);
beta2 = beta2(noZero);
lD = zeros(npar,1);
uD = convexD;
A_eq = [ones(1,npar); y'];
b_eq = [200; 0];
alpha = 1;

%% Optimization Formula %%
% sum(bw*(rho-(alpha*R)./bw).^2) = lambda^2*(alpha*R)^2
% min sum([rho(k)-lambda1/beta1(k)](+)^2*beta2(k)), lambda1=alpha*R, fval=(R-alpha*R)^2
    function [fval, grad] = fn(par)
        p = E' * par;
        M = abs(p) - lambda1./beta1;
        M(M < 0) = 0;
        fval = beta2' * M.^2;
        grad = 2*E * (M.*sign(p).*beta2);
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
        H = E*diag(beta2)*E';
%         options = optimset('Algorithm','interior-point-convex','Display','off');
%         [u, R, ~] = quadprog(H+H', [], [], [], A_eq, b_eq, lD, uD, [], options);
        Prob = qpAssign(H+H', zeros(npar,1), A_eq, b_eq, b_eq, lD, uD, [], 'QP');
        Result = tomRun('qpopt', Prob, 0);
        u = Result.x_k;
        R = Result.f_k;
        x0 = u;
    end;
    
    if(lambda1 > 0)
%         options = optimset('Algorithm','sqp','Display','off','GradObj','on','UseParallel','always');
%         [u, R, ~] = fmincon(@fn, x0, [], [], A_eq, b_eq, lD, uD, [], options);
%         [u1, R1, ~] = snopt(x0,lD,uD,[0;b_eq],[Inf;b_eq],'fnn',[ones(npar,1);y],[ones(npar,1)*2;ones(npar,1)*3],[1:npar,1:npar]',ones(npar,1),(1:npar)');
        Prob = conAssign(@zf, @zg, [], [], lD, uD, 'NLP', x0, [], [], A_eq, b_eq, b_eq);
        Result = tomRun('npsol', Prob, 0);
        u = Result.x_k;
        R = Result.f_k;
    end
    R = lambda1 + sqrt(R/nfea);
    alpha = lambda1/R;
end

%% Solve linear system for lambda1=max (alpha=1) and generating initial point if not provided
if(alpha > 0.999)
    [model, score, alpha] = unifsR_lin(TRAIN, y, TEST, beta1_ori, convexD);
    return;
end

%% Calculate Support Vector, separation hyperplane
tol = 1e-8;
D = max(u)/100;
SV = zeros(1,npar);
SV(u < tol) = -1;
SV(u > uD-tol) = 1;

p = E' * u;
rw = abs(p) - R*alpha./beta1;
rw(rw < 0) = 0;
rw = rw.*sign(p).*beta2;  % sum(rw.^2./beta2/N) - lambda2^2;
rm = rw/(rw'*p)*R;

A_mat = full(E(y==1,:) * rm);
B_mat = full(E(y~=1,:) * rm);
A_sv = SV(y==1);
B_sv = SV(y~=1);
% Calculate opt_a for Group A
if(any(~A_sv))
    opt_a = mean(A_mat(A_sv==0));
else
    opt_a = (min(A_mat(A_sv==1)) + max(A_mat(A_sv==-1)))/2;
end
% Calculate opt_b for Group B
if(any(~B_sv))
    opt_b = -mean(B_mat(B_sv==0));
else
    opt_b = (min(-B_mat(B_sv==-1)) + max(-B_mat(B_sv==1)))/2;
end
dsplane = 100 * (opt_a-opt_b);
rm = [-opt_a-opt_b;2*rm];
w = zeros(nfea+1,1);
w([true,noZero]) = rm/dsplane*100;
delta = 2*R/(opt_a-opt_b)/(rw'*p);

%% Output model and scores %%
model.u = u;
model.R = R;
model.w = w;
model.d = dsplane;
model.D = [D,max(uD)];
model.delta = delta;
model.alpha = alpha;
score.train = [ones(npar,1),TRAIN]*w + 10^-20;
score.test = score.train;
if(~isempty(TEST))
    score.test = [ones(size(TEST,1),1),TEST]*w + 10^-20;
end
end