pj = parpool('local', 16);
addpath(genpath(pwd));
import bioma.data.*
addpath(genpath('/share/apps/MATHWORKS/R2014a/tomlab/'));
pctRunOnAll addpath(genpath('/share/apps/MATHWORKS/R2014a/tomlab/'));
run /share/apps/MATHWORKS/R2014a/tomlab/startup

%% HMP OTU Level Cross Validation Test %%
clearvars -except pj
load('./0.hmp.input/HMPv13.phylo.mat');
load('./0.hmp.input/HMPv13.c19.cvp.mat');

parlist.D = 100./linspace(120,20,30);

parlist.D = 100./linspace(100,20,30);
for k = 1:51
cv = cvp{k};
train = struct();
test = struct();
train.X = HMPv13.edgematrix.mat(cv.train_label,cv.subnode);
train.y = cv.train_group;
test.X = HMPv13.edgematrix.mat(cv.test_label,cv.subnode);
test.y = cv.test_group;
branch = get(HMPv13.phy_tree, 'Distance');
branch = branch(cv.subnode);
foldid = cv.foldid;

tic;
[CV_model, CV_fs, CV_fs1, acctab, score, foldid, time0, CV_alphalist, CV_Dlist, unifs_cv] = unifsR_cv_v1(train, test, foldid, branch, parlist);
time = toc;
disp(time);
save(sprintf('./0.hmp.c19.unifsR0.C/HMPv13.c19.unifsR_v1.%02d.mat',k), 'CV_model', 'CV_fs', 'CV_fs1', 'foldid', 'acctab', 'score', 'time0', 'time', 'CV_alphalist', 'CV_Dlist', 'unifs_cv', '-v7.3');

tic;
[CV_model, CV_fs, CV_fs1, acctab, score, foldid, time0, CV_alphalist, CV_Dlist, unifs_cv] = unifsR_cv_v2(train, test, foldid, branch, parlist);
time = toc;
disp(time);
save(sprintf('./0.hmp.c19.unifsR0.C/HMPv13.c19.unifsR_v2.%02d.mat',k), 'CV_model', 'CV_fs', 'CV_fs1', 'foldid', 'acctab', 'score', 'time0', 'time', 'CV_alphalist', 'CV_Dlist', 'unifs_cv', '-v7.3');

fprintf('c19.%02d finished', k);
clearvars -except i k pj
end


%% Feature Selection Matrix %%
load('./0.hmp.input/HMPv13.phylo.mat');
load('./0.hmp.input/HMPv13.c19.cvp.mat');
for k = 1:51
load(sprintf('./0.hmp.c19.unifsR0.C/HMPv13.c19.unifsR_v1.%02d.mat',k));
cv = cvp{k};
subsample = unique(sort([cv.train_label;cv.test_label]));
subnode = cv.subnode(abs(CV_fs.w(2:end)) > 0);

fsmatrix = DataMatrix(HMPv13.edgebool.mat(:,subnode), HMPv13.edgebool.rownames, HMPv13.edgebool.colnames(subnode));
[~,J] = sort(sum(fsmatrix));
[~,I] = sort(fsmatrix.RowNames);
fsmatrix = fsmatrix(I,J);
for p = size(fsmatrix,2):-1:1
    [~,I] = sort(double(fsmatrix(:,p)), 'descend');
    fsmatrix = fsmatrix(I,:);
end
fsmatrix = DataMatrix(full([CV_fs.w(abs(CV_fs.w)>0)';0,sum(fsmatrix);zeros(HMPv13.Ntaxa,1),double(fsmatrix)]), ['0Coefficient';'0Counts';fsmatrix.RowNames], ['Intercept',fsmatrix.ColNames]);
dmwrite(fsmatrix, sprintf('./0.hmp.c19.unifsR0.C.table/HMPv13.c19.unifsR_v1.%02d.fs.table.%.4f.csv', [k,CV_fs.alpha]), 'Delimiter',',');

%% Distance Matrix %%
N = length(train.subtaxa);
edge_matrix.mat = HMPv13.edgematrix.mat(subsample, subnode);
edge_matrix.rownames = HMPv13.edgematrix.rownames(subsample);
edge_matrix.colnames = HMPv13.edgematrix.colnames(subnode);
edgedist = train.branch(abs(CV_fs.w(2:end)) > 0);
alpha = CV_fs.alpha;
R = CV_fs.d;
alpha0 = 1/(sqrt(N)*(1/alpha-1)+1);
R0 = alpha*R+sqrt(N)*(1-alpha)*R;
dmwrite(DataMatrix(full(edge_matrix.mat), edge_matrix.rownames, edge_matrix.colnames), sprintf('./0.hmp.c19.unifsR0.C.table/HMPv13.c%02d.unifsR0.%02d.core_mat.%.4f.csv', [i,k,CV_fs.alpha]), 'Delimiter', ',');

D0 = pdist(edge_matrix.mat, @(x,J) double((repmat(x,size(J,1),1)-J)).^2*edgedist);
dPCoA = DataMatrix(squareform(sqrt(D0)), edge_matrix.rownames, edge_matrix.rownames);
dmwrite(dPCoA, sprintf('./0.hmp.c19.unifsR0.C.table/HMPv13.c%02d.unifsR0.%02d.fs.dist.dPCoA.%.4f.csv',[i,k,alpha0]), 'Delimiter',',');
D1 = pdist(edge_matrix.mat, @(x,J) abs(double(repmat(x,size(J,1),1)-J))*edgedist);
unifrac = DataMatrix(squareform(D1), edge_matrix.rownames, edge_matrix.rownames);
dmwrite(unifrac, sprintf('./0.hmp.c19.unifsR0.C.table/HMPv13.c%02d.unifsR0.%02d.fs.dist.w.non_nor.%.4f.csv',[i,k,alpha0]), 'Delimiter',',');
Dalpha = pdist(edge_matrix.mat, @(x,J) arrayfun(@(y) unifsR_dis(x,J(y,:),edgedist,edgedist,alpha,N), 1:size(J,1), 'UniformOutput', true));
mydist = DataMatrix(squareform(Dalpha), edge_matrix.rownames, edge_matrix.rownames);
dmwrite(mydist, sprintf('./0.hmp.c19.unifsR0.C.table/HMPv13.c%02d.unifsR0.%02d.fs.dist.alpha.%.4f.csv',[i,k,alpha0]), 'Delimiter',',');
Dfs = alpha0/R0*D1 + (1-alpha0)/R0^2*D0;
mydist = DataMatrix(squareform(Dfs), edge_matrix.rownames, edge_matrix.rownames);
dmwrite(mydist, sprintf('./0.hmp.c19.unifsR0.C.table/HMPv13.c%02d.unifsR0.%02d.fs.dist.comb.%.4f.csv',[i,k,alpha0]), 'Delimiter',',');

clearvars -except i k;
end

delete(pj);

