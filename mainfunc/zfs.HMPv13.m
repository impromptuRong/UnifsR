pj = parpool('local', 16);
addpath(genpath(pwd));
import bioma.data.*

load('./0.hmp.input/HMPv13.phylo.mat');
seqdep = sum(HMPv13.otu_table, 2);
data = HMPv13.otu_table./repmat(seqdep, 1, HMPv13.Ntaxa);
clearvars HMPv13 seqdep

for k = 1:50
disp('***********************************');
%% HMP EdgeMatrix: Classification/FeatureSelection %%
load(sprintf('./0.hmp.input/HMPv13.c19.%02d.mat', k));
% load('./0.hmp.output/HMPv13.c19.unifsR0.00.mat');
fstrain = struct('X', train.edgematrix, 'y', train.bodysite);
fstest = struct('X', test.edgematrix, 'y', test.bodysite);
train.branch(train.branch<1e-4) = 0;
penalty = 1./train.branch;
[edge_acc, edge_feature, edge_model, edge_par, foldid] = cvfs(fstrain, fstest, foldid, penalty, []);

%% HMP OTU Matrix: Classification/FeatureSelection %%
subtaxa = find(mean(data(train.subsample,:), 1)>=0.0001);
fstrain = struct('X', data(train.subsample, subtaxa), 'y', train.bodysite);
fstest = struct('X', data(test.subsample, subtaxa), 'y', test.bodysite);
[otu_acc1, otu_feature1, otu_model1, otu_par1, foldid] = cvfs(fstrain, fstest, foldid, [], []);

subtaxa = find(mean(data(train.subsample,:), 1)>=0.001);
fstrain = struct('X', data(train.subsample, subtaxa), 'y', train.bodysite);
fstest = struct('X', data(test.subsample, subtaxa), 'y', test.bodysite);
[otu_acc2, otu_feature2, otu_model2, otu_par2, foldid] = cvfs(fstrain, fstest, foldid, [], []);

%% Save fs Data and Write Index into Text file %%
save(sprintf('./0.hmp.c19.fs.matlab/HMPv13.c19.fs.R0.%02d.mat', k), 'edge_acc', 'edge_feature', 'edge_model', 'edge_par', 'otu_acc1', 'otu_feature1', 'otu_model1', 'otu_par1', 'otu_acc2', 'otu_feature2', 'otu_model2', 'otu_par2', 'foldid')

fileID = fopen(sprintf('./0.hmp.c19.list/HMPv13.c19.list.R0.%02d.xls', k), 'w');
fprintf(fileID, '%s', 'edgetaxa');
fprintf(fileID, '\t%d', train.subtaxa);
fprintf(fileID, '\n%s', 'otutaxa');
fprintf(fileID, '\t%d', subtaxa);
fprintf(fileID, '\n%s', 'train');
fprintf(fileID, '\t%d', train.subsample);
fprintf(fileID, '\n%s', 'test');
fprintf(fileID, '\t%d', test.subsample);
fprintf(fileID, '\n%s', 'foldid');
fprintf(fileID, '\t%d', foldid);
fclose(fileID)

disp(sprintf('Group.c19.%02d finished', k));
disp('***********************************');
clearvars -except k data pj
end

%% Save Result into Mat file and Txt file%%
for k = 1:50
load(sprintf('./0.hmp.c19.fs.matlab/HMPv13.c19.fs.R0.%02d.mat', k));
colnames = {'No_fea','CVacc','Train_acc','Train_sen','Train_spe','Test_acc','Test_sen','Test_spe','All_acc','All_sen','All_spe'};
rownames = {'svm_linear'; 'svmRFE_linear'; 'svm_rbf'; 'svmRFE_rbf'; 'Logistic'; 'L1reg_Logistic'; 'L1reg_L2loss_SVM'; 'Glmnet_1se'; 'Glmnet_min'; 'Glmnet_1se_w'; 'Glmnet_min_w'};
OutMat = ['edgeMatrix',colnames,'otuMatrix',colnames; rownames,num2cell(edge_acc),rownames,num2cell(otu_acc1)];
fileID = fopen(sprintf('./0.hmp.c19.acc/HMPv13.c19.%02d.acc.xls', k), 'w');
fprintf(fileID, '%s', OutMat{1});
fprintf(fileID, '\t%s', OutMat{1,2:end});
formatSpec = '%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f';
for m = 1:length(rownames)
    fprintf(fileID, ['\n', formatSpec,'\t',formatSpec], OutMat{m+1,:});
end
fclose(fileID);
end
