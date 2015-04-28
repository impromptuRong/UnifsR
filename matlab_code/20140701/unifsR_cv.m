function [CV_model, CV_score, CV_fs, CV_D, CV_alpha, foldid, time0, CVacc, CVauc, TCacc] = unifsR_cv(edgematrix, group, branch, parlist)
%% Turning parameters %%
edgematrix = double(edgematrix);
if isempty(branch)
    branch = ones(size(edgematrix,2),1);
end
beta = branch;
group = nominal(group);
[y, label] = grp2idx(group);
y = 2*y-3;

if(~isfield(parlist,'D') || isempty(parlist.D))
    kD = min([sum(y>0),sum(y<0),18]);
    parlist.D = 1./[linspace(1,sum(y<0)/1.2,kD/1.2);linspace(1,sum(y>0)/1.2,kD/1.2)];
end
if(~isfield(parlist,'kfold') || isempty(parlist.kfold))
    parlist.kfold = 10;
end

if length(parlist.kfold) == length(group)
    foldid = parlist.kfold;
else
    if isnumeric(parlist.kfold)
        cp = cvpartition(group, 'k', parlist.kfold);
    else
        cp = parlist.kfold;
    end
    foldid = sum(cell2mat(arrayfun(@(x) cp.test(x)*x, 1:cp.NumTestSets, 'UniformOutput',false)),2);
end
kfold = max(foldid);
kD = size(parlist.D,2);
kalpha = 101;

%% Generating cross-validation matrix %%
unifs_cv = cell(kD, kfold);
for x = 1:kD
    for k = 1:kfold
	adjD = parlist.D(:,x);
	adjD = parlist.D(:,x)*length(foldid)/sum(foldid~=k);
        unifs_cv{x,k} = {edgematrix(foldid~=k,:);y(foldid~=k);edgematrix(foldid==k,:);y(foldid==k);adjD};
    end
    unifs_cv{x,k+1} = {edgematrix;y;edgematrix;y;parlist.D(:,x)};
end
unifs_cv = reshape(unifs_cv, [], 1);

%% Parallel cross-validation Test %%
tic;
parfor x = 1:size(unifs_cv,1)
fprintf('fold %2d:\tD = %6.2f, %6.2f\n', [ceil(x/kD);1./unifs_cv{x}{5}]);
    model_list = cell(101,1);
    score_list = cell(101,1);
%     alpha_list = zeros(101,1);
%     R_list = zeros(101,1);
    
    % Initial Solver with lambda2=0(alpha=1)
    [model, cv_dec, ~] = unifsR_lin(unifs_cv{x}{1}, unifs_cv{x}{2}, unifs_cv{x}{3}, beta, unifs_cv{x}{5});
    model.label = label;
    score_list(end) = {cv_dec};
    model_list(end) = {model};
%     alpha_list(end) = alpha = 1;
%     R_list(end) = model.R;
    R1 = model.R;
    
    % Initial Solver with lambda1=0(alpha=0)
    [model, cv_dec, ~] = unifsR_con(unifs_cv{x}{1}, unifs_cv{x}{2}, unifs_cv{x}{3}, beta, beta, 0, unifs_cv{x}{5}, model.u);
    model.label = label;
    score_list(1) = {cv_dec};
    model_list(1) = {model};
%     alpha_list(1) = alpha = 0;
%     R_list(1) = model.R;
    R0 = model.R;
    
%    f_alpha = @(x, s) x./(x+s*R0/R1*(1-x));
    f_x = @(alpha, s) 1./(1+R1/R0/s*(1./alpha-1));
    s = 1;    
    for k = 2:100
        Rx = f_x((k-1)/100, s);
        [model, cv_dec, alpha] = unifsR_con(unifs_cv{x}{1}, unifs_cv{x}{2}, unifs_cv{x}{3}, beta, beta, R1*Rx, unifs_cv{x}{5}, model.u);
        model.label = label;
        score_list(k) = {cv_dec};
        model_list(k) = {model};
%         alpha_list(k) = alpha;
%         R_list(k) = model.R;
	s = s * (1/alpha-1)/(100/(k-1)-1);
    end
    unifs_cv{x}{6} = model_list;
    unifs_cv{x}{7} = score_list;
%     unifs_cv{x}{8} = alpha_list;
%     unifs_cv{x}{9} = R_list;
end
time0 = toc;
unifs_cv = reshape(unifs_cv, [], kfold+1);

CV_model = cellfun(@(x) x{6}, unifs_cv(:,end), 'UniformOutput', false);
CV_model = cat(2, CV_model{:});
CV_score = cell(kalpha, kD);
alpha = cellfun(@(x) x.alpha, CV_model, 'UniformOutput', true);
dsmax = cellfun(@(x) x.d, CV_model, 'UniformOutput', true);
% Rmax = cellfun(@(x) x.R, CV_model, 'UniformOutput', true);
Nfea = cellfun(@(x) sum(abs(x.w(2:end))>1e-4), CV_model, 'UniformOutput', true);

CVauc = zeros(kalpha, kD);
CVacc = zeros(kalpha, kD);
TCacc = zeros(kalpha, kD);
parfor x = 1:kD
    cross = unifs_cv(x,1:end-1);
    testy = cell2mat(cellfun(@(x) x{4}', cross, 'UniformOutput', false));
    for k = 1:101
        CV_score{k,x}.tot = [unifs_cv{x,end}{4}, unifs_cv{x,end}{7}{k}.test];
        cv_dec = cell2mat(cellfun(@(x) x{7}{k}.test', cross, 'UniformOutput', false));
        CV_score{k,x}.cv = [testy; cv_dec]';
        [~, ~, ~, auc] = perfcurve(testy, cv_dec, 1);
        CVauc(k,x) = auc;
        confmat = confusionmat(testy, sign(cv_dec));
        CVacc(k,x) = (confmat(1)+confmat(4))/sum(confmat(:));
        confmat = confusionmat(CV_score{k,x}.tot(:,1), sign(CV_score{k,x}.tot(:,2)));
        TCacc(k,x) = (confmat(1)+confmat(4))/sum(confmat(:));
    end
end

score = CVacc*0.8 + TCacc*0.2 - Nfea/size(edgematrix,2)*0.1;
des = sort(CVacc(:),'descend');
cv_index = CVacc >= des(kD+1);
cv_index = find(score >= max(score(cv_index)) & cv_index);
[~,L] = max(mod(cv_index, kalpha));
CV_fs = CV_model{cv_index(L)};
CV_fs.index = cv_index(L);

CV_D = cell(1,kD);
for x = 1:kD
    cv_acc = CVacc(:,x);
    cv_auc = CVauc(:,x);
    tc_acc = TCacc(:,x);
    index = cv_acc >= max(cv_acc)-std(cv_acc)/sqrt(kalpha);
%    index = index & CVauc >= max(CVauc(index))-std(CVauc)/sqrt(kalpha);
    index = index & tc_acc >= max(tc_acc(index));
    index = index & alpha(:,x) >= max(alpha(index, x));
    cvopt.model = CV_model{index, x};
    cvopt.score = CV_score{index, x};
    cvopt.eval = [cv_auc(index), cv_acc(index), tc_acc(index)];
    CV_D(x) = {cvopt};
end

CV_alpha = cell(size(CV_model,1),1);
for x = 1:kalpha
    cv_acc = CVacc(x,:);
    cv_auc = CVauc(x,:);
    tc_acc = TCacc(x,:);
    index = cv_acc>=max(cv_acc)-std(cv_acc)/sqrt(kalpha);
%    index = index & CVauc >= max(CVauc(index))-std(CVauc)/sqrt(kD);
    index = index & tc_acc >= max(tc_acc(index));
    index = index & dsmax(x,:) >= max(dsmax(x,index));
    cvopt.model = CV_model{x, index};
    cvopt.score = CV_score{x, index};
    cvopt.eval = [cv_auc(index), cv_acc(index), tc_acc(index)];
    CV_alpha(x) = {cvopt};
end
end
