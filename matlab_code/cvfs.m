function [matrix_acc, feature_list, model_list, par_list, foldid] = cvfs(train, test, cp, penalty, method)
if ~iscell(train)
    trainX = sparse(double(train.X));
    trainy = 2*grp2idx(train.y) - 3;

    if length(cp) == length(trainy)
        foldid = cp;
    else 
        if isnumeric(cp)
            cp = cvpartition(trainy, 'k', cp);
        end
        foldid = sum(cell2mat(arrayfun(@(x) cp.test(x)*x, 1:cp.NumTestSets, 'UniformOutput',false)),2);
    end
    kfold = max(foldid);
else
    cv_fold = train;
    trainX = sparse(double(cv_fold{end}{1}));
    trainy = 2*grp2idx(cv_fold{end}{2}) - 3;
    kfold = length(cv_fold) - 1;
    foldid = zeros(size(cv_fold{end}{2}));
    for k = 1:kfold
        [~,~,b] = intersect(cv_fold{k}{3}(:,end), cv_fold{end}{3}(:,end));
        foldid(b) = k;
    end
end
if(isempty(method))
    method = {'svm_lin', 'svm_rbf', 'L2reg_Logistic', 'L1reg_Logistic', 'L1reg_L2loss_SVM', 'Glmnet', 'Glmnet_phy'};
end
if(isempty(penalty))
    penalty = ones(size(trainX,2));
end

cv_fold = cell(1, kfold+1);
for k = 1:kfold
    cv_fold{k} = {trainX(foldid~=k,:);trainy(foldid~=k);trainX(foldid==k,:);trainy(foldid==k)};
end
cv_fold{end} = {trainX;trainy;trainX;trainy};
if(~isempty(test) && (isfield(test,'X') && ~isempty(test.X)) && (isfield(test,'y') && ~isempty(test.y)))
    testX = sparse(double(test.X));
    testy = 2*grp2idx(test.y)-3;
else
    testX = trainX;
    testy = trainy;
end

%% Output Result Table %%
matrix_acc = zeros(11, 11);
feature_list = cell(11, 1);
model_list = cell(11, 1);
par_list = cell(11, 1);

%% linear kernel SVM + svmRFE: cost=2^(-5:15) %
if(ismember('svm_lin', method))
    fprintf('\nSVM + svmRFE + linear kernel\n');
    try
        [acc_tab, feature, model] = svmRFE_fold(cv_fold, '-s 0 -v 10 -t 0', [], 'fast');
        % Original Model %
        mod = model{1};
        confmatC = confusionmat(trainy, libsvmpredict(trainy, trainX, mod, '-q'));
        confmatT = confusionmat(testy, libsvmpredict(testy, testX, mod, '-q'));
        confmatA = confmatC + confmatT;
        matrix_acc(1,:) = [acc_tab(1,1:2), [(confmatC(1)+confmatC(4))/sum(confmatC(:)), confmatC(1)/(confmatC(1)+confmatC(3)), confmatC(4)/(confmatC(4)+confmatC(2)), ...
            (confmatT(1)+confmatT(4))/sum(confmatT(:)), confmatT(1)/(confmatT(1)+confmatT(3)), confmatT(4)/(confmatT(4)+confmatT(2)), ...
            (confmatA(1)+confmatA(4))/sum(confmatA(:)), confmatA(1)/(confmatA(1)+confmatA(3)), confmatA(4)/(confmatA(4)+confmatA(2))]*100];
        feature_list(1) = feature(1);
        model_list(1) = {mod};
        par_list(1) = {acc_tab(1,4:6)};        
        % Feature Model %
        opt = find(acc_tab(:,2)==max(acc_tab(:,2)));
        opt = opt(end);
        mod = model{opt};
        confmatC = confusionmat(trainy, libsvmpredict(trainy, trainX(:,feature{opt}), mod, '-q'));
        confmatT = confusionmat(testy, libsvmpredict(testy, testX(:,feature{opt}), mod, '-q'));
        confmatA = confmatC + confmatT;
        matrix_acc(2,:) = [acc_tab(opt,1:2), [(confmatC(1)+confmatC(4))/sum(confmatC(:)), confmatC(1)/(confmatC(1)+confmatC(3)), confmatC(4)/(confmatC(4)+confmatC(2)), ...
            (confmatT(1)+confmatT(4))/sum(confmatT(:)), confmatT(1)/(confmatT(1)+confmatT(3)), confmatT(4)/(confmatT(4)+confmatT(2)), ...
            (confmatA(1)+confmatA(4))/sum(confmatA(:)), confmatA(1)/(confmatA(1)+confmatA(3)), confmatA(4)/(confmatA(4)+confmatA(2))]*100];
        feature_list(2) = feature(opt);
        model_list(2) = {mod};
        par_list(2) = {acc_tab(opt,4:6)};
    catch
        fprintf('SVM + svmRFE + linear kernel Fail !\n');
    end
end

%% RBF kernel SVM + svmRFE: cost=2^(-5:15), gamma=2^(-15:3) %%
if(ismember('svm_rbf', method))
    fprintf('\nSVM + svmRFE + RBF kernel\n');
    try
        [acc_tab, feature, model] = svmRFE_fold(cv_fold, '-s 0 -v 10 -t 2', [], 'fast');
        % Original Model %
        mod = model{1};
        confmatC = confusionmat(trainy, libsvmpredict(trainy, trainX, mod, '-q'));
        confmatT = confusionmat(testy, libsvmpredict(testy, testX, mod, '-q'));
        confmatA = confmatC + confmatT;
        matrix_acc(3,:) = [acc_tab(1,1:2), [(confmatC(1)+confmatC(4))/sum(confmatC(:)), confmatC(1)/(confmatC(1)+confmatC(3)), confmatC(4)/(confmatC(4)+confmatC(2)), ...
            (confmatT(1)+confmatT(4))/sum(confmatT(:)), confmatT(1)/(confmatT(1)+confmatT(3)), confmatT(4)/(confmatT(4)+confmatT(2)), ...
            (confmatA(1)+confmatA(4))/sum(confmatA(:)), confmatA(1)/(confmatA(1)+confmatA(3)), confmatA(4)/(confmatA(4)+confmatA(2))]*100];
        feature_list(3) = feature(1);
        model_list(3) = {mod};
        par_list(3) = {acc_tab(1,4:6)};
        % Feature Model %
        opt = find(acc_tab(:,2)==max(acc_tab(:,2)));
        opt = opt(end);
        mod = model{opt};
        confmatC = confusionmat(trainy, libsvmpredict(trainy, trainX(:,feature{opt}), mod, '-q'));
        confmatT = confusionmat(testy, libsvmpredict(testy, testX(:,feature{opt}), mod, '-q'));
        confmatA = confmatC + confmatT;
        matrix_acc(4,:) = [acc_tab(opt,1:2), [(confmatC(1)+confmatC(4))/sum(confmatC(:)), confmatC(1)/(confmatC(1)+confmatC(3)), confmatC(4)/(confmatC(4)+confmatC(2)), ...
            (confmatT(1)+confmatT(4))/sum(confmatT(:)), confmatT(1)/(confmatT(1)+confmatT(3)), confmatT(4)/(confmatT(4)+confmatT(2)), ...
            (confmatA(1)+confmatA(4))/sum(confmatA(:)), confmatA(1)/(confmatA(1)+confmatA(3)), confmatA(4)/(confmatA(4)+confmatA(2))]*100];
        feature_list(4) = feature(opt);
        model_list(4) = {mod};
        par_list(4) = {acc_tab(opt,4:6)};
    catch
        fprintf('SVM + svmRFE + RBF kernel Fail !\n');
    end
end

%% liblinear L2 logistic Regression: cost=2^(-5:15) %%
if(ismember('L2reg_Logistic', method))
    fprintf('\nL2 Logistic Regression\n');
    try
        log2c = -5:15;
        kfold = 10;
        parfor k = 1:length(log2c)
            cmd = ['-s 0 -v 10 -c ', num2str(2^log2c(k)), ' -q'];
            acck = zeros(kfold, 1);
            for m = 1:kfold
                modk = liblineartrain(cv_fold{m}{2}, cv_fold{m}{1}, regexprep(cmd, '-v (\d+)', ''));
                labelk = liblinearpredict(cv_fold{m}{4}, cv_fold{m}{3}, modk, '-q');
                confmat = confusionmat(cv_fold{m}{4}, labelk);
                acck(m) = (confmat(1)+confmat(4))/sum(confmat(:)) * 100;
            end
            fprintf('Cross Validation Accuracy: %.2f%\n', mean(acck));
            cv(k) = mean(acck);
            md{k} = liblineartrain(cv_fold{end}{2}, cv_fold{end}{1}, regexprep(cmd, '-v (\d+)', ''));
        end
        [acc0, index] = max(cv);
        mod = md{index};
        confmatC = confusionmat(trainy, liblinearpredict(trainy, trainX, mod, '-q'));
        confmatT = confusionmat(testy, liblinearpredict(testy, testX, mod, '-q'));
        confmatA = confmatC + confmatT;
        matrix_acc(5,:) = [size(cv_fold{end}{1},2), acc0, [(confmatC(1)+confmatC(4))/sum(confmatC(:)), confmatC(1)/(confmatC(1)+confmatC(3)), confmatC(4)/(confmatC(4)+confmatC(2)), ...
            (confmatT(1)+confmatT(4))/sum(confmatT(:)), confmatT(1)/(confmatT(1)+confmatT(3)), confmatT(4)/(confmatT(4)+confmatT(2)), ...
            (confmatA(1)+confmatA(4))/sum(confmatA(:)), confmatA(1)/(confmatA(1)+confmatA(3)), confmatA(4)/(confmatA(4)+confmatA(2))]*100];
        feature_list(5) = {1:size(cv_fold{end}{1},2)};
        model_list(5) = {mod};
        par_list(5) = {2^log2c(index)};
        clearvars 'cv';
    catch
        fprintf('L2 Logistic Regression Fail !');
    end
end

%% liblinear L1 Reg. logistic Regression: cost=2^(-5:15) %%
if(ismember('L1reg_Logistic', method))
    fprintf('\nliblinear L1 Reg. logistic Regression\n');
    try
        log2c = -5:15;
        parfor k = 1:length(log2c)
            cmd = ['-s 6 -v 10 -c ', num2str(2^log2c(k)), ' -q'];
            acck = zeros(kfold, 1);
            for m = 1:kfold
                modk = liblineartrain(cv_fold{m}{2}, cv_fold{m}{1}, regexprep(cmd, '-v (\d+)', ''));
                labelk = liblinearpredict(cv_fold{m}{4}, cv_fold{m}{3}, modk, '-q');
                confmat = confusionmat(cv_fold{m}{4}, labelk);
                acck(m) = (confmat(1)+confmat(4))/sum(confmat(:)) * 100;
            end
            fprintf('Cross Validation Accuracy: %.2f%\n', mean(acck));
            cv(k) = mean(acck);
            md{k} = liblineartrain(cv_fold{end}{2}, cv_fold{end}{1}, regexprep(cmd, '-v (\d+)', ''));
        end
        nfea = arrayfun(@(x) sum(abs(md{x}.w)>0), 1:length(log2c), 'UniformOutput', true);
        index = cv==max(cv) & nfea <= min(nfea(cv==max(cv)));
        mod = md{index};
        confmatC = confusionmat(trainy, liblinearpredict(trainy, trainX, mod, '-q'));
        confmatT = confusionmat(testy, liblinearpredict(testy, testX, mod, '-q'));
        confmatA = confmatC + confmatT;
        matrix_acc(6,:) = [sum(abs(mod.w)>0), max(cv), [(confmatC(1)+confmatC(4))/sum(confmatC(:)), confmatC(1)/(confmatC(1)+confmatC(3)), confmatC(4)/(confmatC(4)+confmatC(2)), ...
            (confmatT(1)+confmatT(4))/sum(confmatT(:)), confmatT(1)/(confmatT(1)+confmatT(3)), confmatT(4)/(confmatT(4)+confmatT(2)), ...
            (confmatA(1)+confmatA(4))/sum(confmatA(:)), confmatA(1)/(confmatA(1)+confmatA(3)), confmatA(4)/(confmatA(4)+confmatA(2))]*100];
        feature_list(6) = {find(abs(mod.w)>0)};
        model_list(6) = {mod};
        par_list(6) = {2^log2c(index)};
        clearvars 'cv' 'md';
    catch
        fprintf('L1 Logistic Regression Fail !');
    end
end

%% liblinear L1 Reg. L2 loss SVM: cost=2^(-5:15) %%
if(ismember('L1reg_L2loss_SVM', method))
    fprintf('\nliblinear L1 Reg. L2 loss SVM\n');
    try
        log2c = -5:15;
        parfor k = 1:length(log2c)
            cmd = ['-s 5 -v 10 -c ', num2str(2^log2c(k)), ' -q'];
            acck = zeros(kfold, 1);
            for m = 1:kfold
                modk = liblineartrain(cv_fold{m}{2}, cv_fold{m}{1}, regexprep(cmd, '-v (\d+)', ''));
                labelk = liblinearpredict(cv_fold{m}{4}, cv_fold{m}{3}, modk, '-q');
                confmat = confusionmat(cv_fold{m}{4}, labelk);
                acck(m) = (confmat(1)+confmat(4))/sum(confmat(:)) * 100;
            end
            fprintf('Cross Validation Accuracy: %.2f%\n', mean(acck));
            cv(k) = mean(acck);
            md{k} = liblineartrain(cv_fold{end}{2}, cv_fold{end}{1}, regexprep(cmd, '-v (\d+)', ''));
        end
        nfea = arrayfun(@(x) sum(abs(md{x}.w)>0), 1:length(log2c), 'UniformOutput', true);
        index = cv==max(cv) & nfea <= min(nfea(cv==max(cv)));
        mod = md{index};
        confmatC = confusionmat(trainy, liblinearpredict(trainy, trainX, mod, '-q'));
        confmatT = confusionmat(testy, liblinearpredict(testy, testX, mod, '-q'));
        confmatA = confmatC + confmatT;
        matrix_acc(7,:) = [sum(abs(mod.w)>0), max(cv), [(confmatC(1)+confmatC(4))/sum(confmatC(:)), confmatC(1)/(confmatC(1)+confmatC(3)), confmatC(4)/(confmatC(4)+confmatC(2)), ...
            (confmatT(1)+confmatT(4))/sum(confmatT(:)), confmatT(1)/(confmatT(1)+confmatT(3)), confmatT(4)/(confmatT(4)+confmatT(2)), ...
            (confmatA(1)+confmatA(4))/sum(confmatA(:)), confmatA(1)/(confmatA(1)+confmatA(3)), confmatA(4)/(confmatA(4)+confmatA(2))]*100];
        feature_list(7) = {find(abs(mod.w)>0)};
        model_list(7) = {mod};
        par_list(7) = {2^log2c(index)};
        clearvars 'cv' 'md';
    catch
        fprintf('liblinear L1 Reg. L2 loss SVM Fail !\n');
    end
end

%% Glmnet Logistic Regression: alpha=0:0.01:1 %%
if(ismember('Glmnet', method))
    fprintf('\nGlmnet Logistic Regression\n');
    try
        alpha = 0:0.01:1;
        parfor k = 1:length(alpha)
            opts = struct;
            opts.alpha = alpha(k);
            opts.standardize = false;
            mod = cvglmnet(trainX, trainy, 'binomial', glmnetSet(opts), 'class', 10, foldid);
            md{k} = mod;
            cv_1se(k) = 1-mod.cvm(mod.lambda==mod.lambda_1se);
            cv_min(k) = 1-mod.cvm(mod.lambda==mod.lambda_min);
        end
        % Use Lambda_1se for selection %
        index = find(cv_1se==max(cv_1se));
        index = index(end);
        acc0 = cv_1se(index);
        mod = md{index};
        coef = cvglmnetCoef(mod, 'lambda_1se');
        confmatC = confusionmat(trainy, sign(cvglmnetPredict(mod, trainX, 'lambda_1se')));
        confmatT = confusionmat(testy, sign(cvglmnetPredict(mod, testX, 'lambda_1se')));
        confmatA = confmatC + confmatT;
        matrix_acc(8,:) = [(sum(abs(coef)>0)-1), [acc0, (confmatC(1)+confmatC(4))/sum(confmatC(:)), confmatC(1)/(confmatC(1)+confmatC(3)), confmatC(4)/(confmatC(4)+confmatC(2)), ...
            (confmatT(1)+confmatT(4))/sum(confmatT(:)), confmatT(1)/(confmatT(1)+confmatT(3)), confmatT(4)/(confmatT(4)+confmatT(2)), ...
            (confmatA(1)+confmatA(4))/sum(confmatA(:)), confmatA(1)/(confmatA(1)+confmatA(3)), confmatA(4)/(confmatA(4)+confmatA(2))]*100];
        feature_list(8) = {find(abs(coef(2:end))>0)};
        model_list(8) = {mod};
        par_list(8) = {[alpha(index),mod.lambda_1se]};
        % Use Lambda_min for selection %
        index = find(cv_min==max(cv_min));
        index = index(end);
        acc0 = cv_min(index);
        mod = md{index};
        coef = cvglmnetCoef(mod, 'lambda_min');
        confmatC = confusionmat(trainy, sign(cvglmnetPredict(mod, trainX, 'lambda_min')));
        confmatT = confusionmat(testy, sign(cvglmnetPredict(mod, testX, 'lambda_min')));
        confmatA = confmatC + confmatT;
        matrix_acc(9,:) = [(sum(abs(coef)>0)-1), [acc0, (confmatC(1)+confmatC(4))/sum(confmatC(:)), confmatC(1)/(confmatC(1)+confmatC(3)), confmatC(4)/(confmatC(4)+confmatC(2)), ...
            (confmatT(1)+confmatT(4))/sum(confmatT(:)), confmatT(1)/(confmatT(1)+confmatT(3)), confmatT(4)/(confmatT(4)+confmatT(2)), ...
            (confmatA(1)+confmatA(4))/sum(confmatA(:)), confmatA(1)/(confmatA(1)+confmatA(3)), confmatA(4)/(confmatA(4)+confmatA(2))]*100];
        feature_list(9) = {find(abs(coef(2:end))>0)};
        model_list(9) = {mod};
        par_list(9) = {[alpha(index),mod.lambda_min]};
        clearvars 'cv_1se' 'cv_min' 'md';
    catch
        fprintf('Glmnet Logistic Regression Fail !\n');
    end
end

%% Glmnet Logistic Regression: alpha=0:0.01:1 with branch length %%
if(ismember('Glmnet_phy', method))
    fprintf('\nGlmnet Logistic Regression with penalty_factor\n');
    try
        alpha = 0:0.1:1;
        exclude = isinf(penalty);
        cvTrain = trainX(:, ~exclude);
        cvPenalty = penalty(~exclude);
        parfor k = 1:length(alpha)
            opts = struct;
            opts.alpha = alpha(k);
            opts.standardize = false;
            opts.penalty_factor = cvPenalty;
            % opts.exclude = find(exclude);
            mod = cvglmnet(cvTrain, trainy, 'binomial', glmnetSet(opts), 'class', 10, foldid);
            md{k} = mod;
            cv_1se(k) = 1-mod.cvm(mod.lambda==mod.lambda_1se);
            cv_min(k) = 1-mod.cvm(mod.lambda==mod.lambda_min);
        end
        % Use Lambda_1se for selection %
        index = find(cv_1se==max(cv_1se));
        index = index(end);
        acc0 = cv_1se(index);
        mod = md{index};
        coef = cvglmnetCoef(mod, 'lambda_1se');
        confmatC = confusionmat(trainy, sign(cvglmnetPredict(mod, trainX(:,~exclude), 'lambda_1se')));
        confmatT = confusionmat(testy, sign(cvglmnetPredict(mod, testX(:,~exclude), 'lambda_1se')));
        confmatA = confmatC + confmatT;
        matrix_acc(10,:) = [(sum(abs(coef)>0)-1), [acc0, (confmatC(1)+confmatC(4))/sum(confmatC(:)), confmatC(1)/(confmatC(1)+confmatC(3)), confmatC(4)/(confmatC(4)+confmatC(2)), ...
            (confmatT(1)+confmatT(4))/sum(confmatT(:)), confmatT(1)/(confmatT(1)+confmatT(3)), confmatT(4)/(confmatT(4)+confmatT(2)), ...
            (confmatA(1)+confmatA(4))/sum(confmatA(:)), confmatA(1)/(confmatA(1)+confmatA(3)), confmatA(4)/(confmatA(4)+confmatA(2))]*100];
        feature_list(10) = {find(abs(coef(2:end))>0)};
        model_list(10) = {mod};
        par_list(10) = {[alpha(index),mod.lambda_1se]};
        % Use Lambda_min for selection %
        index = find(cv_min==max(cv_min));
        index = index(end);
        acc0 = cv_min(index);
        mod = md{index};
        coef = cvglmnetCoef(mod, 'lambda_min');
        confmatC = confusionmat(trainy, sign(cvglmnetPredict(mod, trainX(:,~exclude), 'lambda_min')));
        confmatT = confusionmat(testy, sign(cvglmnetPredict(mod, testX(:,~exclude), 'lambda_min')));
        confmatA = confmatC + confmatT;
        matrix_acc(11,:) = [(sum(abs(coef)>0)-1), [acc0, (confmatC(1)+confmatC(4))/sum(confmatC(:)), confmatC(1)/(confmatC(1)+confmatC(3)), confmatC(4)/(confmatC(4)+confmatC(2)), ...
            (confmatT(1)+confmatT(4))/sum(confmatT(:)), confmatT(1)/(confmatT(1)+confmatT(3)), confmatT(4)/(confmatT(4)+confmatT(2)), ...
            (confmatA(1)+confmatA(4))/sum(confmatA(:)), confmatA(1)/(confmatA(1)+confmatA(3)), confmatA(4)/(confmatA(4)+confmatA(2))]*100];
        feature_list(11) = {find(abs(coef(2:end))>0)};
        model_list(11) = {mod};
        par_list(11) = {[alpha(index),mod.lambda_min]};
        clearvars 'cv_1se' 'cv_min' 'md';
    catch
        fprintf('Glmnet Logistic Regression with penalty_factor Fail !\n');
    end
end
end


%% Random Forest %%
% fprintf('\nRandom Forest\n');
% B = TreeBagger(NTrees,X,Y)
% [labels,score] = predict(B,X)
% [acc_tab, feature, model] = rfRFE_fold(cv_fold, 500, [], [], 'fast');
% Original Model %
% mod = model{1};
% pred = classRF_predict(testX, mod);
% confmat = confusionmat(testy, sign(pred));
% acc = (confmat(1)+confmat(4))/sum(confmat(:));
% matrix_acc(5,:) = [acc_tab(1,1:3),acc(1)];
% feature_list(5) = feature(1);
% model_list(5) = {mod};
% par_list(5) = {acc_tab(1,4:6)};
% Feature Model %
% opt = find(acc_tab(:,2)==max(acc_tab(:,2)));
% opt = opt(end);
% mod = model{opt};
% pred = classRF_predict(testX, mod);
% confmat = confusionmat(testy, sign(pred));
% acc = (confmat(1)+confmat(4))/sum(confmat(:));
% matrix_acc(6,:) = [acc_tab(opt,1:3),acc(1)];
% feature_list(6) = feature(opt);
% model_list(6) = {mod};
% par_list(6) = {acc_tab(opt,4:6)};