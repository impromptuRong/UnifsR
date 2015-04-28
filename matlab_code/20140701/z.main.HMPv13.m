matlabpool local 12;
addpath(genpath(pwd));
import bioma.data.*
% addpath(genpath('/storage/scratch/rr0311'));
% pctRunOnAll addpath(genpath('/storage/scratch/rr0311/'));
addpath(genpath('/share/apps/MATHWORKS/R2013b/tomlab/'));
pctRunOnAll addpath(genpath('/share/apps/MATHWORKS/R2013b/tomlab/'));
run /share/apps/MATHWORKS/R2013b/tomlab/startup
% run /home/rr0311/tomlab/startup


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
    metaInfo.(regexprep(field, {'("|'')','[^a-zA-Z0-9]'}, {'','_'})) = tmp;
end
clearvars tmp i field metainfo;
metaInfo.Sam_ID = cellfun(@num2str, num2cell(metaInfo.SampleID), 'UniformOutput', false);

%% Read in Taxonomy Info %%
taxainfo = csvimport('./0.raw/otu_taxa_psn_v13.csv', 'delimiter', ',');
taxaInfo = struct;
for i=1:size(taxainfo,2)
    colinfo = regexprep(taxainfo(:, i), {'("|'')','[^a-zA-Z0-9]'}, {'','_'});
    taxaInfo.(char(colinfo(1))) = colinfo(2:end);
end
clearvars i colinfo taxainfo;

%% Read in Data Info %%
spa = spconvert(load('./0.raw/otu_table_psn_v13.sp.mt'));
dataInfo = full(spa);
rownames = csvimport('./0.raw/otu_table_psn_v13.sp.rn', 'noHeader', true);
colnames = csvimport('./0.raw/otu_table_psn_v13.sp.cn', 'noHeader', true);
colnames = cellfun(@num2str, colnames, 'UniformOutput', false);
dataInfo = DataMatrix(dataInfo, 'RowNames', regexprep(rownames,'[^a-zA-Z0-9]','_'), 'ColNames', regexprep(colnames,'[^a-zA-Z0-9]','_'));
dataInfo = dataInfo';
clearvars spa rownames colnames;

%% Combine Info and Calculate Unifrac Structure %%
HMPv13 = phyloseq('otu_table', dataInfo, 'sam_table', metaInfo, 'phy_tree', './0.raw/rep_set_v13.tre', 'tax_table', taxaInfo);
save('./0.hmp.input/HMPv13.phylo.mat', 'HMPv13', '-v7.3');

%% Generate Select Groups for Validation %%
load('./0.hmp.input/HMPv13.phylo.mat');
% Delete sample depth smaller than 500
phylo = HMPv13;
subsam = sum(phylo.otu_table,2)>=500;
bodysite = phylo.sam_table.HMPbodysubsite;
branchDis = get(phylo.phy_tree, 'Distance');
%% Input Grouping Fiels %%
readinfo = csvimport('HMP.group.csv', 'delimiter', ',');
groupinfo = struct;
for i=1:size(readinfo,2)
    field = regexprep(char(readinfo(1,i)), {'("|'')', '[^a-zA-Z0-9]'}, {'','_'});
    tmp = readinfo(2:end,i);
    if ~iscellstr(tmp)
        tmp = cell2mat(tmp);
    else
        tmp = regexprep(tmp, '("|'')', '');
    end
    groupinfo.(field) = tmp;
end

for k = 1:10
% for i = 1:size(groupinfo.Group1,1)
for i = [19,20,28,32,33,35,39]
    % Find required Group-pairs. Remove low-abundance Taxa. %
    subindex = subsam.*(cellfun(@isempty, regexpi(bodysite, groupinfo.Group1{i})) - cellfun(@isempty, regexpi(bodysite, groupinfo.Group2{i})));
    subtaxa = find(mean(phylo.edgematrix.mat(subindex~=0,:),1) >= 0.001);
%     subtaxa = mean(phylo.otu_table(subindex~=0,:)./repmat(sum(phylo.otu_table(subindex~=0,:),2),1,phylo.Ntaxa), 1)>=0.01;
    % Training Data normalization %
    subtrain = [randsample(find(subindex==-1), groupinfo.G1_v13_train(i)); randsample(find(subindex==1), groupinfo.G2_v13_train(i))];
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
    save(sprintf('./0.hmp.input/HMPv13.c%02d.%02d.mat', [i,k]), 'train', 'test', '-v7.3');
end

end

%% HMP OTU Level Cross Validation Test %%
clear;
for i = 2:45
load(sprintf('./0.hmp.input/HMPv13.c%02d.mat',i));
tic;
[CV_model, CV_score, CV_opt, time0] = unifsR_cv(train.edgematrix*100, train.bodysite, train.branch, []);
time = toc;
disp(time);
[CV_pred, CV_test] = cellfun(@(x) unifsR_pred(x, test.edgematrix*100, 'cell'), CV_model, 'UniformOutput', false);
[opt_pred, opt_test] = cellfun(@(x) unifsR_pred(x.model, test.edgematrix*100, 'cell'), CV_opt, 'UniformOutput', false);
save(sprintf('./0.hmp.input/HMPv13.c%02d.unifsR.mat',i), 'CV_model', 'CV_score', 'CV_opt', 'CV_pred', 'CV_test', 'opt_pred', 'opt_test', 'time0', 'time');
end



load('./0.hmp.input/HMPv13.c19.00.mat');
load('./0.hmp.output/HMPv13.c19.unifsR1.00.mat');
train0 = train;
test0 = test;

[~, score] = unifsR_pred(CV_fs, train0.edgematrix, 'nominal');
decidata = [grp2idx(test0.bodysite)*2-3, score, 1./(1+exp(-score))-0.5];
bpart1 = decidata(:,1)==1 & decidata(:,2)<1;
gpart1 = decidata(:,1)==1 & decidata(:,2)>=1;
bpart2 = decidata(:,1)~=1 & decidata(:,2)>-1;
gpart2 = decidata(:,1)~=1 & decidata(:,2)<=-1;
% [f1, x1] = ksdensity(decidata(decidata(:,1)~=1,3), 'bandwidth', 0.02);
% [f2, x2] = ksdensity(decidata(decidata(:,1)==1,3), 'bandwidth', 0.02);
% figure(1)
% hold on;
% plot(x1, f1, 'Col', 'red');
% plot(x2, f2, 'Col', 'blue');

i = 19;
for k = 1:50
[subtrain,idx] = sort([randsample(find(bpart1),80); randsample(find(gpart1),40); randsample(find(bpart2),80); randsample(find(gpart2),40)]);
subtest = setdiff(1:size(decidata,1), subtrain);
foldid = [mod(randsample(1:80,80),10),mod(randsample(1:40,40),10),mod(randsample(1:80,80),10),mod(randsample(1:40,40),10)]+1;
foldid = foldid(idx);

train = struct();
train.edgematrix = train0.edgematrix(subtrain,:);
train.bodysite = train0.bodysite(subtrain,:);
train.branch = train0.branch;
train.subsample = train0.subsample(subtrain);
train.subtaxa = train0.subtaxa;

test = struct();
test.edgematrix = train0.edgematrix(subtest,:);
test.bodysite = train0.bodysite(subtest,:);
test.branch = train0.branch;
test.subsample = train0.subsample(subtest);
test.subtaxa = train0.subtaxa;
save(sprintf('./0.hmp.input/HMPv13.c%02d.%02d.mat', [i,k]), 'train', 'test', 'foldid');
end


