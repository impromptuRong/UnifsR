matlabpool local 12;
addpath(genpath(pwd));
import bioma.data.*
addpath(genpath('/share/apps/MATHWORKS/R2013b/tomlab/'));
pctRunOnAll addpath(genpath('/share/apps/MATHWORKS/R2013b/tomlab/'));
run /share/apps/MATHWORKS/R2013b/tomlab/startup

phylo = phyloseq('otu_table', './0.raw/Skin/skin.summary.otu.csv', 'phy_tree', './0.raw/Skin/skin.otu.tre.trans', 'sam_table', './0.raw/Skin/skin.site.meta.csv', 'tax_table', './0.raw/Skin/skin.dic.otu.csv', 'seqdep', true);
save('./0.matlab.out/skin_otu.phylo.mat', 'phylo');
phytreewrite('./0.matlab.out/skin_otu.ref.tre', phylo.phy_tree);
% dmwrite(DataMatrix(full(phylo.edgebool.mat), phylo.edgebool.rownames, phylo.edgebool.colnames), './0.matlab.out/skin_otu.ref.tre.edge.csv', 'Delimiter', ',');
dmwrite(DataMatrix(full(phylo.edgematrix.mat), phylo.edgematrix.rownames, phylo.edgematrix.colnames), './0.matlab.out/skin_otu.edgematrix.csv', 'Delimiter', ',');

% Delete sample depth smaller than 500
bodysite = phylo.sam_table.Sample_site;
branchDis = get(phylo.phy_tree, 'Distance');
%% Input Grouping Fiels %%
groupinfo.Group1 = {'palm', 'palm', 'plantar_foot'};
groupinfo.Group2 = {'plantar_foot', 'popliteal_fossa', 'popliteal_fossa'};
groupinfo.train1 = [40, 40, 40];
groupinfo.train2 = [40, 30, 30];

for i = 1:size(groupinfo.Group1,1)
    % Find required Group-pairs. Remove low-abundance Taxa. %
    subindex = cellfun(@isempty, regexpi(bodysite, groupinfo.Group1{i})) - cellfun(@isempty, regexpi(bodysite, groupinfo.Group2{i}));
    subtaxa = find(mean(phylo.edgematrix.mat(subindex~=0,:),1) >= 0.001);
    % Training Data normalization %
    subtrain = [randsample(find(subindex==-1), groupinfo.train1(i)); randsample(find(subindex==1), groupinfo.train2(i))];
    train = struct();
    train.edgematrix = phylo.edgematrix.mat(subtrain, subtaxa);
    train.bodysite = nominal(subindex(subtrain), {groupinfo.Group1{i},groupinfo.Group2{i}});
    train.branch = branchDis(subtaxa);
    train.subsample = subtrain;
    train.subtaxa = subtaxa;
    
    subtest = setdiff(find(subindex~=0), subtrain);
    test = struct();
    test.edgematrix = phylo.edgematrix.mat(subtest, subtaxa);
    test.bodysite = nominal(subindex(subtest), {groupinfo.Group1{i},groupinfo.Group2{i}});
    test.branch = branchDis(subtaxa);
    test.subsample = subtest;
    test.subtaxa = subtaxa;
    save(sprintf('./0.hmp.input/skin.c%02d.mat', i), 'train', 'test');
end

%% HMP OTU Level Cross Validation Test %%
clear;
i = 1;
load(sprintf('./0.hmp.input/skin.c%02d.mat',i));
parlist.kfold=10;
parlist.D = 1./[1:2:30; 1:2:30];
tic;
[CV_model, CV_score, CV_fs, CV_D, CV_alpha, cp, time0, CVacc, CVauc, TCacc] = unifsR_cv(train.edgematrix, train.bodysite, train.branch, []);
time = toc;
disp(time);
[CV_pred, CV_test] = cellfun(@(x) unifsR_pred(x, test.edgematrix, 'cell'), CV_model, 'UniformOutput', false);
[D_pred, D_test] = cellfun(@(x) unifsR_pred(x.model, test.edgematrix, 'cell'), CV_D, 'UniformOutput', false);
[alpha_pred, alpha_test] = cellfun(@(x) unifsR_pred(x.model, test.edgematrix, 'cell'), CV_alpha, 'UniformOutput', false);
TPacc = cellfun(@(x) confusionmat(test.bodysite,nominal(x)), CV_pred, 'UniformOutput', false);
TPacc = cellfun(@(x) (x(1)+x(4))/sum(x(:)), TPacc, 'UniformOutput', true);
[CVacc(:,2),TCacc(:,2),TPacc(:,2)]

save(sprintf('./0.hmp.output/skin.c%02d.unifsR.mat',i), 'CV_model', 'CV_score', 'CV_fs', 'CV_D', 'CV_alpha', 'CVacc', 'CVauc', 'TCacc', 'TPacc', 'CV_pred', 'CV_test', 'D_pred', 'D_test', 'alpha_pred', 'alpha_test', 'time0', 'time', 'cp');

