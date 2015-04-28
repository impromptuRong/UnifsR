function [CV_model, CV_score, CV_fs, foldid, time0, CVacc, CVauc, TCacc, CV_alphalist, CV_Dlist, unifs_cv] = unifsR_cv(edgematrix, group, branch, parlist)
% pj = parpool('local', 16);
% addpath(genpath(pwd));
% import bioma.data.*
% addpath(genpath('/share/apps/MATHWORKS/R2014a/tomlab/'));
% pctRunOnAll addpath(genpath('/share/apps/MATHWORKS/R2014a/tomlab/'));
% run /share/apps/MATHWORKS/R2014a/tomlab/startup
% clearvars -except pj;
% load('./0.matlab.out/oral_otu.phylo.mat');
% edgematrix = train.edgematrix;
% group = train.group;
% branch = train.branch;
% parlist=[];
% i=19;
% k=0;
% load(sprintf('./0.hmp.input/HMPv13.c%02d.%02d.mat',[i,k]));
% edgematrix = train.edgematrix*100;
% group = train.bodysite;
% branch = train.branch;
% parlist.kfold = foldid;
% parlist.D = 100./linspace(120,20,30);
% if k~=0
%     parlist.D = 100./linspace(90,20,30);
% end


%% Checking Input Data %%
edgematrix = double(edgematrix);
if isempty(branch)
    branch = ones(size(edgematrix,2),1);
end
beta = branch;
group = nominal(group);
[y, label] = grp2idx(group);
y = 2*y-3;
kalpha = 101;

%% Analyzing Cross-Validation ID %%
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

%% Generate Loss Parameters Seq %%
if(~isfield(parlist,'D') || isempty(parlist.D))
    model0 = unifsR_con(edgematrix, y, [], beta, beta, 0, 100*abs(y), []);
    model1 = unifsR_lin(edgematrix, y, [], beta, 100*abs(y));
    parlist.D = 100./linspace(floor((length(y)-abs(sum(y)))*0.4*(1-1/kfold)), floor(200/max([model0.u;model1.u])), 20);
%     parlist.C = 2.^((0:20)-round(log2(quantile(max(edgematrix,[],2),0.99))));
end
kD = length(parlist.D);
Dlist = parlist.D;

%% Generate models for different parameters %%
unifs_all = cell(kD, 1);
parfor x = 1:kD
fprintf('fold 00:\tD = %8.3f\n', Dlist(x));
    convexA = Dlist(x) * abs(y);    
    model_list = cell(kalpha,1);
    score_list = cell(kalpha,1);
    alpha_list = zeros(kalpha,1);
    R_list = zeros(kalpha,1);
    D_list = zeros(kalpha,1);

    %  Calculate model1.C for alpha = 1 (lambda2 = max)
    [model, cv_dec, ~] = unifsR_lin(edgematrix, y, edgematrix, beta, convexA);
    model.label = label;
    score_list(end) = {cv_dec};
    model_list(end) = {model};
    alpha_list(end) = model.alpha;
    R_list(end) = model.R;
    D_list(end) = model.D(1);
    R1 = model.R;
    
    %  Calculate model0.C for alpha = 0 (lambda2 = 0)
    [model, cv_dec, ~] = unifsR_con(edgematrix, y, edgematrix, beta, beta, 0, convexA, []);
    model.label = label;
    score_list(1) = {cv_dec};
    model_list(1) = {model};
    alpha_list(1) = model.alpha;
    R_list(1) = model.R;
    D_list(1) = model.D(1);
    R0 = model.R;
 
    % Use R0 and R1 to generate lambda sequence
    f_x = @(alpha, s) 1./(1+R1/R0/s*(1./alpha-1));  % f_alpha = @(x, s) x./(x+s*R0/R1*(1-x));
    s = 1;
    model = model_list{1};
    plotdata = zeros(kalpha,5);
    plotdata(1,:) = [0,1,0,0,R0];
    plotdata(end,:) = [1,1,R1,1,R1];
    for k = 2:kalpha-1
        Rx = f_x((k-1)/(kalpha-1), s);
        [model, cv_dec, ~] = unifsR_con(edgematrix, y, edgematrix, beta, beta, Rx*R1, convexA, model.u);
        model.label = label;
        score_list(k) = {cv_dec};
        model_list(k) = {model};
        
        plotdata(k,:) = [Rx,s,R1*Rx,model.alpha,model.R];
        alpha_list(k) = model.alpha;
        R_list(k) = model.R;
        D_list(k) = model.D(1);
	    s = s * (1/model.alpha-1)/((kalpha-1)/(k-1)-1);
    end
%     save('tmplot.mat','plotdata');

    unifs_all{x} = {edgematrix; y; edgematrix; y};
    unifs_all{x}{5} = [model_list{1}.D(2)*model_list{1}.delta, model_list{kalpha}.D(2)*model_list{kalpha}.delta, Dlist(x)];
    
    unifs_all{x}{6} = model_list;
    unifs_all{x}{7} = score_list;
    unifs_all{x}{8} = alpha_list;
    unifs_all{x}{9} = R_list;
    unifs_all{x}{10} = D_list;
end
modelC = cellfun(@(x) [x{6}{1}.D(2)*x{6}{1}.delta, x{6}{kalpha}.D(2)*x{6}{kalpha}.delta, x{6}{1}.D(2)], unifs_all, 'UniformOutput', false);

%% Generating cross-validation Dataset %%
unifs_cv = cell(kD, kfold);
for x = 1:kD
    for k = 1:kfold
        unifs_cv{x,k} = {edgematrix(foldid~=k,:);y(foldid~=k);edgematrix(foldid==k,:);y(foldid==k);[modelC{x}(1:2),modelC{x}(3)]};
        % /sum(foldid~=k)*length(y)
    end
end
unifs_cv = reshape(unifs_cv', [], 1);
% [cell2mat(cellfun(@(x) x{5}, unifs_cv,'UniformOutput',false)),cell2mat(modelC(td(:),:))]

%% Parallel cross-validation Test for models %%
tic;
% cvDlist = cell(size(unifs_cv));
parfor x = 1:size(unifs_cv,1)
fprintf('fold %2d:\tD = %8.3f\n', [rem(x, kfold); unifs_cv{x}{5}(3)]);
    model_list = cell(kalpha,1);
    score_list = cell(kalpha,1);
    alpha_list = zeros(kalpha,1);
    R_list = zeros(kalpha,1);
    D_list = zeros(kalpha,1);

    % Generate convex-hull constrain: convexD1 based on model.C1
    [model, cv_dec, convexD1] = unifsR_init1(unifs_cv{x}{1}, unifs_cv{x}{2}, unifs_cv{x}{3}, beta, unifs_cv{x}{5}(2));
    model.label = label;
    score_list(end) = {cv_dec};
    model_list(end) = {model};
    alpha_list(end) = model.alpha;
    R_list(end) = model.R;
    D_list(end) = model.D(1);
    R1 = model.R;
    
    %  Generate convex-hull constrain: convexD0 based on model.C0
    [model, cv_dec, convexD0] = unifsR_init0(unifs_cv{x}{1}, unifs_cv{x}{2}, unifs_cv{x}{3}, beta, unifs_cv{x}{5}(1));
    model.label = label;
    score_list(1) = {cv_dec};
    model_list(1) = {model};
    alpha_list(1) = model.alpha;
    R_list(1) = model.R;
    D_list(1) = model.D(1);
    R0 = model.R;
    
%     cvDlist(x) = {[max(convexD0), max(convexD1)]};
% end
% cvDlist = [reshape(cvDlist, [], kfold), num2cell(Dlist')];
% save('tmplot.mat','cvDlist');

    % Use R0 and R1 to generate lambda sequence
    f_x = @(alpha, s) 1./(1+R1/R0/s*(1./alpha-1));  % f_alpha = @(x, s) x./(x+s*R0/R1*(1-x));
    s = 1;
    model = model_list{1};
    plotdata = zeros(kalpha,5);
    plotdata(1,:) = [0,1,0,0,R0];
    plotdata(end,:) = [1,1,R1,1,R1];
    for k = 2:kalpha-1
        Rx = f_x((k-1)/(kalpha-1), s);
        convexD = convexD1*(k-1)/100 + convexD0*(1-(k-1)/100);
        [model, cv_dec, ~] = unifsR_con(unifs_cv{x}{1}, unifs_cv{x}{2}, unifs_cv{x}{3}, beta, beta, Rx*R1, convexD, model.u);
        model.label = label;
        score_list(k) = {cv_dec};
        model_list(k) = {model};
        
        plotdata(k,:) = [Rx,s,R1*Rx,model.alpha,model.R];
        alpha_list(k) = model.alpha;
        R_list(k) = model.R;
        D_list(k) = model.D(1);
	    s = s * (1/model.alpha-1)/((kalpha-1)/(k-1)-1);
    end
%     save('tmplot.mat','plotdata');
    
    unifs_cv{x}{6} = model_list;
    unifs_cv{x}{7} = score_list;
    unifs_cv{x}{8} = alpha_list;
    unifs_cv{x}{9} = R_list;
    unifs_cv{x}{10} = D_list;
end
time0 = toc;
unifs_cv = [reshape(unifs_cv, kfold, [])', unifs_all];

CV_model = cellfun(@(x) x{6}, unifs_cv(:,end), 'UniformOutput', false);
CV_model = cat(2, CV_model{:});
CV_score = cell(kalpha, kD);
% alpha = cellfun(@(x) x.alpha, CV_model, 'UniformOutput', true);
% dsmax = cellfun(@(x) x.d, CV_model, 'UniformOutput', true);
% Rmax = cellfun(@(x) x.R, CV_model, 'UniformOutput', true);
Nfea = cellfun(@(x) sum(abs(x.w(2:end))>1e-4), CV_model, 'UniformOutput', true);
CV_alphalist = cellfun(@(x) x{8}, unifs_cv, 'UniformOutput', false);
CV_Dlist = cellfun(@(x) x{10}, unifs_cv, 'UniformOutput', false);

CVauc = zeros(kalpha, kD);
CVacc = zeros(kalpha, kD);
TCacc = zeros(kalpha, kD);
for x = 1:kD
    cross = unifs_cv(x,1:end-1);
    testy = cell2mat(cellfun(@(x) x{4}', cross, 'UniformOutput', false));
    for k = 1:kalpha
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

% CV_D = cell(1,kD);
% for x = 1:kD
%     cv_acc = CVacc(:,x);
%     cv_auc = CVauc(:,x);
%     tc_acc = TCacc(:,x);
%     index = cv_acc >= max(cv_acc)-std(cv_acc)/sqrt(kalpha);
% %    index = index & CVauc >= max(CVauc(index))-std(CVauc)/sqrt(kalpha);
%     index = index & tc_acc >= max(tc_acc(index));
%     index = index & alpha(:,x) >= max(alpha(index, x));
%     index = find(index,1,'last');
%     cvopt.model = CV_model{index, x};
%     cvopt.score = CV_score{index, x};
%     cvopt.eval = [cv_auc(index), cv_acc(index), tc_acc(index)];
%     CV_D(x) = {cvopt};
% end
% 
% CV_alpha = cell(size(CV_model,1),1);
% for x = 1:kalpha
%     cv_acc = CVacc(x,:);
%     cv_auc = CVauc(x,:);
%     tc_acc = TCacc(x,:);
%     index = cv_acc>=max(cv_acc)-std(cv_acc)/sqrt(kalpha);
% %    index = index & CVauc >= max(CVauc(index))-std(CVauc)/sqrt(kD);
%     index = index & tc_acc >= max(tc_acc(index));
%     index = index & dsmax(x,:) >= max(dsmax(x,index));
%     index = find(index,1,'last');
%     cvopt.model = CV_model{x, index};
%     cvopt.score = CV_score{x, index};
%     cvopt.eval = [cv_auc(index), cv_acc(index), tc_acc(index)];
%     CV_alpha(x) = {cvopt};
% end
end
