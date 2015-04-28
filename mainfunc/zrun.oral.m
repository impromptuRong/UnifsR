pj = parpool('local', 16);
addpath(genpath(pwd));
import bioma.data.*
addpath(genpath('/share/apps/MATHWORKS/R2014a/tomlab/'));
pctRunOnAll addpath(genpath('/share/apps/MATHWORKS/R2014a/tomlab/'));
run /share/apps/MATHWORKS/R2014a/tomlab/startup

%% Oral Level Cross Validation Test %%
otu = phyloseq('otu_table', './0.raw/Oral/oral.summary.otu.unrounded.csv', 'sam_table', './0.raw/Oral/meta.info.csv', 'phy_tree', './0.raw/Oral/oral.otu.tre', 'tax_table', './0.raw/Oral/dic.otu.csv');
genus = phyloseq('otu_table', './0.raw/Oral/oral.summary.genus.unrounded.csv', 'sam_table', './0.raw/Oral/meta.info.csv', 'phy_tree', './0.raw/Oral/oral.genus.tre.trans', 'tax_table', './0.raw/Oral/dic.genus.csv');
phylum = phyloseq('otu_table', './0.raw/Oral/oral.summary.phylum.unrounded.csv', 'sam_table', './0.raw/Oral/meta.info.csv', 'phy_tree', './0.raw/Oral/oral.phylum.tre.trans', 'tax_table', './0.raw/Oral/dic.phylum.csv');
oralnew = phyloseq('otu_table', './0.raw/Oral/oralnew.summary.otu.unrounded.csv', 'sam_table', './0.raw/Oral/oralnew.meta.info.csv', 'phy_tree', './0.raw/Oral/oralnew.otu.tre', 'tax_table', './0.raw/Oral/oralnew.dic.otu.csv');

parlist.size = [1;1];
cvp = cv_partition(phylum, 'Periodontitis', {'P-','P+'}, parlist);
foldid1 = cvp.foldid;
cvp = cv_partition(oralnew, 'Group', [], parlist);
foldid2 = cvp.foldid;
save('./oral.phylo.mat', 'otu','genus','phylum','oralnew','foldid1','foldid2');

%% Oral Phylum %%
load('./0.matlab.out/oral.phylo.mat');
phylo = phylum;
train = struct();
train.X = phylo.edgematrix.mat;
train.y = phylo.sam_table.Periodontitis;
branch = get(phylo.phy_tree, 'Distance');

tic;
[CV_model, CV_fs, CV_fs1, acctab, score, foldid, time0, CV_alphalist, CV_Dlist, unifs_cv] = unifsR_cv_v1(train, [], foldid1, branch, []);
time = toc;
disp(time);
save('./0.matlab.out/oral_phylum.unifsR_v1.mat', 'CV_model', 'CV_fs', 'CV_fs1', 'foldid', 'acctab', 'score', 'time0', 'time', 'CV_alphalist', 'CV_Dlist');

%% Oral Genus %%
load('./0.matlab.out/oral.phylo.mat');
phylo = genus;
train = struct();
train.X = phylo.edgematrix.mat;
train.y = phylo.sam_table.Periodontitis;
branch = get(phylo.phy_tree, 'Distance');

tic;
[CV_model, CV_fs, CV_fs1, acctab, score, foldid, time0, CV_alphalist, CV_Dlist, unifs_cv] = unifsR_cv_v1(train, [], foldid1, branch, []);
time = toc;
disp(time);
save('./0.matlab.out/oral_genus.unifsR_v1.mat', 'CV_model', 'CV_fs', 'CV_fs1', 'foldid', 'acctab', 'score', 'time0', 'time', 'CV_alphalist', 'CV_Dlist');

%% Oral otu %%
load('./0.matlab.out/oral.phylo.mat');
phylo = otu;
train = struct();
train.X = phylo.edgematrix.mat;
train.y = phylo.sam_table.Periodontitis;
branch = get(phylo.phy_tree, 'Distance');

tic;
[CV_model, CV_fs, CV_fs1, acctab, score, foldid, time0, CV_alphalist, CV_Dlist, unifs_cv] = unifsR_cv_v1(train, [], foldid1, branch, []);
time = toc;
disp(time);
save('./0.matlab.out/oral_otu.unifsR_v1.mat', 'CV_model', 'CV_fs', 'CV_fs1', 'foldid', 'acctab', 'score', 'time0', 'time', 'CV_alphalist', 'CV_Dlist');

%% Oral new %%
load('./0.matlab.out/oral.phylo.mat');
phylo = oralnew;
train = struct();
train.X = phylo.edgematrix.mat;
train.y = phylo.sam_table.Group;
branch = get(phylo.phy_tree, 'Distance');

tic;
[CV_model, CV_fs, CV_fs1, acctab, score, foldid, time0, CV_alphalist, CV_Dlist, unifs_cv] = unifsR_cv_v1(train, [], foldid2, branch, []);
time = toc;
disp(time);
save('./0.matlab.out/oral_new.unifsR_v1.mat', 'CV_model', 'CV_fs', 'CV_fs1', 'foldid', 'acctab', 'score', 'time0', 'time', 'CV_alphalist', 'CV_Dlist');

delete(gcp);


