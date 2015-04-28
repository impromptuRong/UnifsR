matlabpool local 8;
addpath(genpath(pwd));
import bioma.data.*
% addpath(genpath('/storage/scratch/rr0311/matlab_code/'));
% pctRunOnAll addpath(genpath('/storage/scratch/rr0311/matlab_code/'));
% poolobj = gcp;
% addAttachedFiles(poolobj,{'/storage/scratch/rr0311/matlab_code/'});

load('./0.hmp.input/HMPv13.phylo.mat')
load('./0.hmp.input/HMPv13.c19.mat')
% load('./0.hmp.output/HMPv13.c19.unifsR.mat')

%% HMP EdgeMatrix %%
edgematrix_acc = zeros(7, 3);

% linear kernel grid search: cost=2^(-5:15) %
log2c = -5:15;
parfor k = 1:length(log2c)
    cmd = ['-s 0 -t 0 -v 10 -c ', num2str(2^log2c(k)), ' -q'];
    cv(k) = libsvmtrain(2*grp2idx(train.bodysite)-3, double(train.edgematrix), cmd);
end
[acc0, index] = max(cv);
% acc0 = libsvmtrain(2*grp2idx(train.bodysite)-3, double(train.edgematrix), ['-s 0 -t 0 -v 10 -c ', num2str(2^log2c(index)), ' -q']);
mod = libsvmtrain(2*grp2idx(train.bodysite)-3, double(train.edgematrix), ['-s 0 -t 0 -c ', num2str(2^log2c(index)), ' -q']);
[~,acc1,~] = libsvmpredict(2*grp2idx(train.bodysite)-3, double(train.edgematrix), mod, '-q');
[~,acc2,~] = libsvmpredict(2*grp2idx(test.bodysite)-3, double(test.edgematrix), mod, '-q');
edgematrix_acc(1,:) = [acc0,acc1(1),acc2(1)];
clearvars 'cv';

% rbf kernel grid search: gamma=2^(-15:3),cost=2^(-5:15) % 
log2c = -5:15;
log2g = -15:3;
[g1, g2] = meshgrid(log2c, log2g);
exp_grid = [g1(:), g2(:)];
parfor k = 1:size(exp_grid,1)
    cmd = ['-s 0 -t 2 -v 10 -c ', num2str(2^exp_grid(k,1)), ' -g ', num2str(2^exp_grid(k,2)), ' -q'];
    cv(k) = libsvmtrain(2*grp2idx(train.bodysite)-3, double(train.edgematrix), cmd);
end
% fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', exp_grid(index,1), exp_grid(index,2), cv(index), bestc, bestg, bestcv);
[acc0, index] = max(cv);
% acc0 = libsvmtrain(2*grp2idx(train.bodysite)-3, double(train.edgematrix), ['-s 0 -t 2 -v 10 -c ', num2str(2^exp_grid(index,1)), ' -g ', num2str(2^exp_grid(index,2)), ' -q']);
mod = libsvmtrain(2*grp2idx(train.bodysite)-3, double(train.edgematrix), ['-s 0 -t 2 -c ', num2str(2^exp_grid(index,1)), ' -g ', num2str(2^exp_grid(index,2)), ' -q']);
[~,acc1,~] = libsvmpredict(2*grp2idx(train.bodysite)-3, double(train.edgematrix), mod, '-q');
[~,acc2,~] = libsvmpredict(2*grp2idx(test.bodysite)-3, double(test.edgematrix), mod, '-q');
edgematrix_acc(2,:) = [acc0,acc1(1),acc2(1)];
clearvars 'cv';

% liblinear L2 logistic Regression: cost=2^(-5:15) %
log2c = -5:15;
parfor k = 1:length(log2c)
    cmd = ['-s 0 -v 10 -c ', num2str(2^log2c(k)), ' -q'];
    cv(k) = liblineartrain(2*grp2idx(train.bodysite)-3, train.edgematrix, cmd);
end
[acc0, index] = max(cv);
% acc0 = liblineartrain(2*grp2idx(train.bodysite)-3, train.edgematrix, ['-s 0 -v 10 -c ', num2str(2^log2c(index)), ' -q']);
mod = liblineartrain(2*grp2idx(train.bodysite)-3, train.edgematrix, ['-s 0 -c ', num2str(2^log2c(index)), ' -q']);
[~,acc1,~] = liblinearpredict(2*grp2idx(train.bodysite)-3, train.edgematrix, mod, '-q');
[~,acc2,~] = liblinearpredict(2*grp2idx(test.bodysite)-3, test.edgematrix, mod, '-q');
edgematrix_acc(3,:) = [acc0,acc1(1),acc2(1)];
clearvars 'cv';

% liblinear L1 Reg. logistic Regression: cost=2^(-5:15) %
log2c = -5:15;
parfor k = 1:length(log2c)
    cmd = ['-s 6 -v 10 -c ', num2str(2^log2c(k)), ' -q'];
    cv(k) = liblineartrain(2*grp2idx(train.bodysite)-3, train.edgematrix, cmd);
end
[acc0, index] = max(cv);
% acc0 = liblineartrain(2*grp2idx(train.bodysite)-3, train.edgematrix, ['-s 6 -v 10 -c ', num2str(2^log2c(index)), ' -q']);
mod = liblineartrain(2*grp2idx(train.bodysite)-3, train.edgematrix, ['-s 6 -c ', num2str(2^log2c(index)), ' -q']);
[~,acc1,~] = liblinearpredict(2*grp2idx(train.bodysite)-3, train.edgematrix, mod, '-q');
[~,acc2,~] = liblinearpredict(2*grp2idx(test.bodysite)-3, test.edgematrix, mod, '-q');
edgematrix_acc(4,:) = [acc0,acc1(1),acc2(1)];
clearvars 'cv';

% liblinear L1 Reg. L2 loss SVM: cost=2^(-5:15) %
log2c = -5:15;
parfor k = 1:length(log2c)
    cmd = ['-s 5 -v 10 -c ', num2str(2^log2c(k)), ' -q'];
    cv(k) = liblineartrain(2*grp2idx(train.bodysite)-3, train.edgematrix, cmd);
end
[acc0, index] = max(cv);
% acc0 = liblineartrain(2*grp2idx(train.bodysite)-3, train.edgematrix, ['-s 5 -v 10 -c ', num2str(2^log2c(index)), ' -q']);
mod = liblineartrain(2*grp2idx(train.bodysite)-3, train.edgematrix, ['-s 5 -c ', num2str(2^log2c(index)), ' -q']);
[~,acc1,~] = liblinearpredict(2*grp2idx(train.bodysite)-3, train.edgematrix, mod, '-q');
[~,acc2,~] = liblinearpredict(2*grp2idx(test.bodysite)-3, test.edgematrix, mod, '-q');
edgematrix_acc(5,:) = [acc0,acc1(1),acc2(1)];
clearvars 'cv';

% Glmnet Logistic Regression: cost=2^(-5:15) %
alpha = 0:0.1:1;
parfor k = 1:length(alpha)
    opts = struct;
    opts.alpha = alpha(k);
    mod = cvglmnet(train.edgematrix, 2*grp2idx(train.bodysite)-3, 'binomial', glmnetSet(opts), 'class', 10);
    enet{k} = mod;
    cv(k) = 1-mod.cvm(mod.lambda==mod.lambda_1se);
end
index = find(cv==max(cv));
acc0 = cv(index);
pred = cvglmnetPredict(enet{index}, train.edgematrix, 'lambda_1se');
confmat = confusionmat(2*grp2idx(train.bodysite)-3, sign(pred));
acc1 = (confmat(1)+confmat(4))/sum(confmat(:));
pred = cvglmnetPredict(enet{index}, test.edgematrix, 'lambda_1se');
confmat = confusionmat(2*grp2idx(test.bodysite)-3, sign(pred));
acc2 = (confmat(1)+confmat(4))/sum(confmat(:));
edgematrix_acc(7,:) = [acc0,acc1,acc2];

% Random Forest %
% B = TreeBagger(NTrees,X,Y)
% [labels,score] = predict(B,X)
mod = classRF_train(train.edgematrix, 2*grp2idx(train.bodysite)-3, 1000, 5);
pred = classRF_predict(train.edgematrix, mod);
confmat = confusionmat(2*grp2idx(train.bodysite)-3, sign(pred));
acc1 = (confmat(1)+confmat(4))/sum(confmat(:));
pred = classRF_predict(test.edgematrix, mod);
confmat = confusionmat(2*grp2idx(test.bodysite)-3, sign(pred));
acc2 = (confmat(1)+confmat(4))/sum(confmat(:));
edgematrix_acc(6,:) = [1-model.errtr(end,1),acc1,acc2];


%% HMP OTU Matrix %%
otutable_acc = zeros(7, 3);
seqdep = sum(HMPv13.otu_table, 2);
data = HMPv13.otu_table./repmat(seqdep, 1, HMPv13.Ntaxa);
nozero = mean(data(train.subsample,:), 1)>=0.001;
traindata = double(data(train.subsample, nozero));
testdata = double(data(test.subsample, nozero));

% linear kernel grid search: cost=2^(-5:15) %
log2c = -5:15;
parfor k = 1:length(log2c)
    cmd = ['-s 0 -t 0 -v 10 -c ', num2str(2^log2c(k)), ' -q'];
    cv(k) = libsvmtrain(2*grp2idx(train.bodysite)-3, traindata, cmd);
end
[acc0, index] = max(cv);
% acc0 = libsvmtrain(2*grp2idx(train.bodysite)-3, traindata, ['-s 0 -t 0 -v 10 -c ', num2str(2^log2c(index)), ' -q']);
mod = libsvmtrain(2*grp2idx(train.bodysite)-3, traindata, ['-s 0 -t 0 -c ', num2str(2^log2c(index)), ' -q']);
[~,acc1,~] = libsvmpredict(2*grp2idx(train.bodysite)-3, traindata, mod, '-q');
[~,acc2,~] = libsvmpredict(2*grp2idx(test.bodysite)-3, testdata, mod, '-q');
otutable_acc(1,:) = [acc0,acc1(1),acc2(1)];
clearvars 'cv';

% rbf kernel grid search: gamma=2^(-15:3),cost=2^(-5:15) % 
log2c = -5:15;
log2g = -15:3;
[g1, g2] = meshgrid(log2c, log2g);
exp_grid = [g1(:), g2(:)];
parfor k = 1:size(exp_grid,1)
    cmd = ['-s 0 -t 2 -v 10 -c ', num2str(2^exp_grid(k,1)), ' -g ', num2str(2^exp_grid(k,2)), ' -q'];
    cv(k) = libsvmtrain(2*grp2idx(train.bodysite)-3, traindata, cmd);
end
% fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', exp_grid(index,1), exp_grid(index,2), cv(index), bestc, bestg, bestcv);
[acc0, index] = max(cv);
% acc0 = libsvmtrain(2*grp2idx(train.bodysite)-3, traindata, ['-s 0 -t 2 -v 10 -c ', num2str(2^exp_grid(index,1)), ' -g ', num2str(2^exp_grid(index,2)), ' -q']);
mod = libsvmtrain(2*grp2idx(train.bodysite)-3, traindata, ['-s 0 -t 2 -c ', num2str(2^exp_grid(index,1)), ' -g ', num2str(2^exp_grid(index,2)), ' -q']);
[~,acc1,~] = libsvmpredict(2*grp2idx(train.bodysite)-3, traindata, mod, '-q');
[~,acc2,~] = libsvmpredict(2*grp2idx(test.bodysite)-3, testdata, mod, '-q');
otutable_acc(2,:) = [acc0,acc1(1),acc2(1)];
clearvars 'cv';

% liblinear L2 logistic Regression: cost=2^(-5:15) %
log2c = -5:15;
parfor k = 1:length(log2c)
    cmd = ['-s 0 -v 10 -c ', num2str(2^log2c(k)), ' -q'];
    cv(k) = liblineartrain(2*grp2idx(train.bodysite)-3, sparse(traindata), cmd);
end
[acc0, index] = max(cv);
% acc0 = liblineartrain(2*grp2idx(train.bodysite)-3, sparse(traindata), ['-s 0 -v 10 -c ', num2str(2^log2c(index)), ' -q']);
mod = liblineartrain(2*grp2idx(train.bodysite)-3, sparse(traindata), ['-s 0 -c ', num2str(2^log2c(index)), ' -q']);
[~,acc1,~] = liblinearpredict(2*grp2idx(train.bodysite)-3, sparse(traindata), mod, '-q');
[~,acc2,~] = liblinearpredict(2*grp2idx(test.bodysite)-3, sparse(testdata), mod, '-q');
otutable_acc(3,:) = [acc0,acc1(1),acc2(1)];
clearvars 'cv';

% liblinear L1 Reg. logistic Regression: cost=2^(-5:15) %
log2c = -5:15;
parfor k = 1:length(log2c)
    cmd = ['-s 6 -v 10 -c ', num2str(2^log2c(k)), ' -q'];
    cv(k) = liblineartrain(2*grp2idx(train.bodysite)-3, sparse(traindata), cmd);
end
[acc0, index] = max(cv);
% acc0 = liblineartrain(2*grp2idx(train.bodysite)-3, sparse(traindata), ['-s 6 -v 10 -c ', num2str(2^log2c(index)), ' -q']);
mod = liblineartrain(2*grp2idx(train.bodysite)-3, sparse(traindata), ['-s 6 -c ', num2str(2^log2c(index)), ' -q']);
[~,acc1,~] = liblinearpredict(2*grp2idx(train.bodysite)-3, sparse(traindata), mod, '-q');
[~,acc2,~] = liblinearpredict(2*grp2idx(test.bodysite)-3, sparse(testdata), mod, '-q');
otutable_acc(4,:) = [acc0,acc1(1),acc2(1)];
clearvars 'cv';

% liblinear L1 Reg. L2 loss SVM: cost=2^(-5:15) %
log2c = -5:15;
parfor k = 1:length(log2c)
    cmd = ['-s 5 -v 10 -c ', num2str(2^log2c(k)), ' -q'];
    cv(k) = liblineartrain(2*grp2idx(train.bodysite)-3, sparse(traindata), cmd);
end
[acc0, index] = max(cv);
% acc0 = liblineartrain(2*grp2idx(train.bodysite)-3, sparse(traindata), ['-s 5 -v 10 -c ', num2str(2^log2c(index)), ' -q']);
mod = liblineartrain(2*grp2idx(train.bodysite)-3, sparse(traindata), ['-s 5 -c ', num2str(2^log2c(index)), ' -q']);
[~,acc1,~] = liblinearpredict(2*grp2idx(train.bodysite)-3, sparse(traindata).edgematrix, mod, '-q');
[~,acc2,~] = liblinearpredict(2*grp2idx(test.bodysite)-3, sparse(testdata).edgematrix, mod, '-q');
otutable_acc(5,:) = [acc0,acc1(1),acc2(1)];
clearvars 'cv';

% Glmnet Logistic Regression: cost=2^(-5:15) %
alpha = 0:0.1:1;
parfor k = 1:length(alpha)
    opts = struct;
    opts.alpha = alpha(k);
    mod = cvglmnet(sparse(traindata), 2*grp2idx(train.bodysite)-3, 'binomial', glmnetSet(opts), 'class', 10);
    enet{k} = mod;
    cv(k) = 1-mod.cvm(mod.lambda==mod.lambda_1se);
end
index = find(cv==max(cv));
acc0 = cv(index);
pred = cvglmnetPredict(enet{index}, sparse(traindata), 'lambda_1se');
confmat = confusionmat(2*grp2idx(train.bodysite)-3, sign(pred));
acc1 = (confmat(1)+confmat(4))/sum(confmat(:));
pred = cvglmnetPredict(enet{index}, sparse(testdata), 'lambda_1se');
confmat = confusionmat(2*grp2idx(test.bodysite)-3, sign(pred));
acc2 = (confmat(1)+confmat(4))/sum(confmat(:));
otutable_acc(7,:) = [acc0,acc1,acc2];

% Random Forest %
% B = TreeBagger(NTrees,X,Y)
% [labels,score] = predict(B,X)
model = classRF_train(sparse(double(HMPv13.otu_table(train.subsample,:))*100), 2*grp2idx(train.bodysite)-3, 500);





