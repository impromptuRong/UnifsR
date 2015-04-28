matlabpool local 12;
addpath(genpath(pwd));
import bioma.data.*

%% HMP EdgeMatrix: Classification/FeatureSelection %%
load('./0.hmp.input/HMPv13.c19.mat');
load('./0.hmp.output/HMPv13.c19.unifsR2.mat');
fstrain = struct('X', train.edgematrix, 'y', train.bodysite);
fstest = struct('X', test.edgematrix, 'y', test.bodysite);
clearvars -except fstrain fstest cp
[edge_acc, edge_feature, edge_model, edge_par, foldid] = cvfs(fstrain, fstest, cp);

%% HMP OTU Matrix: Classification/FeatureSelection %%
load('./0.hmp.input/HMPv13.c19.mat');
load('./0.hmp.input/HMPv13.phylo.mat');
seqdep = sum(HMPv13.otu_table, 2);
data = HMPv13.otu_table./repmat(seqdep, 1, HMPv13.Ntaxa);
subtaxa = find(mean(data(train.subsample,:), 1)>=0.001);
fstrain = struct('X', data(train.subsample, subtaxa), 'y', train.bodysite);
fstest = struct('X', data(test.subsample, subtaxa), 'y', test.bodysite);
[otu_acc, otu_feature, otu_model, otu_par, foldid] = cvfs(fstrain, fstest, foldid);

%% Save Result into Mat file %%
save('./0.hmp.output/HMPv13.c19.fs.mat', 'edge_acc', 'edge_feature', 'edge_model', 'edge_par', 'otu_acc', 'otu_feature', 'otu_model', 'otu_par', 'foldid')
%% Write Index into Text file %%
fileID = fopen('./0.hmp.output/HMPv13.c19.list.txt', 'w');
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
fclose(fileID);


