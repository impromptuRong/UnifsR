matlabpool local 8;
addpath(genpath(pwd));
import bioma.data.*

%% Read in Meta Info %%
metainfo = csvimport('./0.raw/v13_map_uniquebyPSN.txt', 'delimiter', '\t');
metaInfo = struct;
for i=1:size(metainfo,2)
    tmp = metainfo(2:end, i);
    if ~iscellstr(tmp)
        tmp = cell2mat(tmp);
    end
    field = char(metainfo(1, i));
    if(isempty(field))
        field = strcat('NA',num2str(i));
    end
    metaInfo.(regexprep(field,'[^a-zA-Z0-9]','_')) = tmp;
end
clearvars tmp i field metainfo;
metaInfo.SampleID = cellfun(@num2str, num2cell(metaInfo.SampleID), 'UniformOutput', false);

%% Read in Taxonomy Info %%
taxainfo = csvimport('./0.raw/otu_taxa_psn_v13.csv', 'delimiter', ',');
taxaInfo = struct;
for i=1:size(taxainfo,2)
    colinfo = regexprep(taxainfo(:, i), '[^a-zA-Z0-9]', '_');
    taxaInfo.(char(colinfo(1))) = colinfo(2:end);
end
clearvars i colinfo taxainfo;

%% Read in Tree Info %%
tree = phytreeread('./0.raw/rep_set_v13.tre');
tree = reroot(tree);
leafname = strrep(get(tree,'LeafNames'), '''', '');
leafname = regexprep(leafname, '[^a-zA-Z0-9]', '_');
branch = get(tree, 'Distance');
branch(branch<=0) = 1e-8;
tree = phytree(get(tree,'Pointers'), branch, leafname);
clearvars leafname branch;

%% Read in Data Info %%
spa = spconvert(load('./0.raw/otu_table_psn_v13.sp.mt'));
dataInfo = full(spa);
rownames = csvimport('./0.raw/otu_table_psn_v13.sp.rn', 'noHeader', true);
colnames = csvimport('./0.raw/otu_table_psn_v13.sp.cn', 'noHeader', true);
colnames = cellfun(@num2str, colnames, 'UniformOutput', false);
dataInfo = DataMatrix(dataInfo, 'RowNames', regexprep(rownames,'[^a-zA-Z0-9]','_'), 'ColNames', regexprep(colnames,'[^a-zA-Z0-9]','_'));
dataInfo = dataInfo';
clearvars spa rownames colnames;

%% Combine Info and Calculate Unifrac Structure%%
phylo = phyloseq('otu_table', dataInfo, 'sam_table', metaInfo, 'phy_tree', tree, 'tax_table', taxaInfo);
save('./0.matlab.out/HMPv1v3.phylo.mat', 'phylo');

%% HMP OTU Level Cross Validation Test %%
% load('./0.matlab.out/HMP_100_Att_Sup.phylo.mat');
% phylo_train = phyloseq('otu_table', sub.train.data, 'phy_tree', sub.tree, 'sam_table', sub.train.meta, 'tax_table', sub.taxa);
% save('./0.matlab.out/HMP_100_Att_Sup.phylo_train.mat', 'phylo_train');
% phylo_test = phyloseq('otu_table', sub.test.data, 'phy_tree', sub.tree, 'sam_table', sub.test.meta, 'tax_table', sub.taxa);
% save('./0.matlab.out/HMP_100_Att_Sup.phylo_test.mat', 'phylo_test');
% phytreewrite('./0.matlab.out/HMP_100_Att_Sup.ref.tre', phylo_train.phy_tree);
% dmwrite(DataMatrix(full(phylo_train.edgebool.mat), phylo_train.edgebool.rownames, phylo_train.edgebool.colnames), './0.matlab.out/HMP_100_Att_Sup.ref.tre.edge.csv', 'Delimiter', ',');
% dmwrite([phylo_train.otu_table; phylo_test.otu_table], './0.matlab.out/HMP_100_Att_Sup.summary.csv', 'Delimiter', ',');
% clearvars sub;

clear;
load('./0.matlab.out/HMP_100_Att_Sup.phylo_train.mat');
load('./0.matlab.out/HMP_100_Att_Sup.phylo_test.mat');
parlist.kfold = 10;
parlist.D = 1./(1:2:11);
tic;
[CV_model, CV_score, CV_opt, time0] = unifsR_cv(phylo_train.edgematrix.mat*100, phylo_train.sam_table.HMPbodysubsite, get(phylo_train.phy_tree, 'Distance'), parlist);
time = toc;
disp(time);
[CV_pred, CV_test] = cellfun(@(x) unifsR_pred(x, phylo_test.edgematrix.mat*100, 'cell'), CV_model, 'UniformOutput', false);
[opt_pred, opt_test] = cellfun(@(x) unifsR_pred(x.model, phylo_test.edgematrix.mat*100, 'cell'), CV_opt, 'UniformOutput', false);
save('./0.matlab.out/HMP_100_Att_Sup.unifsR.mat', 'CV_model', 'CV_score', 'CV_opt', 'CV_pred', 'CV_test', 'opt_pred', 'opt_test', 'time0', 'time');
confusionmat(CV_opt{6}.score.cv(:,1),sign(CV_opt{6}.score.cv(:,2)))
[rX,rY,~,auc] = perfcurve(CV_opt{6}.score.cv(:,1),CV_opt{6}.score.cv(:,2),1);
confusionmat(phylo_test.sam_table.HMPbodysubsite, opt_pred{6});

% libSVM Analysis %
X = double(phylo_train.otu_table);
edgematrix = full(phylo_train.edgematrix.mat);
group = phylo_train.sam_table.HMPbodysubsite;
gclass = class(group);
[y, label] = grp2idx(group);
y = 2*y-3;

SVMmodel = libsvmtrain(y, X, '-s 0 -t 0');
SVMmodel_cv = libsvmtrain(y, X, '-s 0 -t 0 -v 10');

plotroc(y, X, SVMmodel);
plotroc(y, X, '-s 0 -t 0 -v 10');
plotroc(y, edgematrix, '-s 0 -t 0 -v 10');
