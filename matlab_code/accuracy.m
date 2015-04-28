addpath(genpath(pwd));
import bioma.data.*

%% Load HMPv13.c19.unifsR0.k.mat %%
list = load('list.csv');
N = length(list);
result = zeros(N, 11);
colnames = {'No_fea','CVacc','Train_acc','Train_sen','Train_spe','Test_acc','Test_sen','Test_spe','All_acc','All_sen','All_spe'};
rownames = {'svm_linear'; 'svmRFE_linear'; 'svm_rbf'; 'svmRFE_rbf'; 'Logistic'; 'L1reg_Logistic'; 'L1reg_L2loss_SVM'; 'Glmnet_1se'; 'Glmnet_min'};

for p = 1:N
    k = list(p);
    load(sprintf('./0.hmp.input/HMPv13.c19.%02d.mat',k));
    load(sprintf('./0.hmp.c19.unifsR0/HMPv13.c19.unifsR0.%02d.mat',k));
    m = CV_fs.index;
    confmatC = confusionmat(CV_score{m}.tot(:,1), sign(CV_score{m}.tot(:,2)));
    confmatT = confusionmat(test.bodysite, nominal(CV_pred{m}));
    confmatA = confmatC + confmatT;
    result(p,:) = [sum(abs(CV_fs.w)>0)-1, [CVacc(m), TCacc(m), confmatC(1)/(confmatC(1)+confmatC(3)), confmatC(4)/(confmatC(4)+confmatC(2)), ...
        TPacc(m), confmatT(1)/(confmatT(1)+confmatT(3)), confmatT(4)/(confmatT(4)+confmatT(2)), ...
        (confmatA(1)+confmatA(4))/sum(confmatA(:)), confmatA(1)/(confmatA(1)+confmatA(3)), confmatA(4)/(confmatA(4)+confmatA(2))]*100];

    OutMat = ['unifsR0',num2cell(result(p,:)),'unifsR0',num2cell(result(p,:))];
    fileID = fopen(sprintf('./0.hmp.c19.acc/HMPv13.c19.%02d.acc.xls', k), 'a');
    formatSpec = '%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f';
    fprintf(fileID, [formatSpec,'\t',formatSpec], OutMat{:});
    fclose(fileID);
end

csvwrite('HMPv13.c19.fs.unifsR0.edgeacc.csv', result);

for p = 1:N
    k = list(p);
    load(sprintf('./0.hmp.output/HMPv13.c19.fs.R0.%02d.mat',k));
    fstab = [uftab*100; edge_acc; otu_acc1; otu_acc2];
    rf = csvread(sprintf('./0.hmp.output/HMPv13.c19.fs.R0.%02d.mat',k));
    fstab([],:) = rf;
    result(k) = fstab;
end








clear;
load('./0.hmp.input/HMPv13.phylo.mat')
phytreewrite('./0.raw/HMP/HMPv13.ref.tre', HMPv13.phy_tree);
mmwrite('./0.raw/HMPv13.otu.mtx', sparse(double(HMPv13.otu_table)));
mmwrite('./0.raw/HMP/HMPv13.edge.mtx', HMPv13.edgematrix.mat);

save('./temp.mat', 'CV_model', 'CV_score', 'CV_fs', 'foldid', 'CVacc', 'CVauc', 'TCacc')

%% Plot model selection figures %%
load('./0.hmp.input/HMPv13.phylo.mat')
load('./0.hmp.output/HMPv13.c19.unifsR2.mat')
[kalpha, kD] = size(CV_model);
alpha = cellfun(@(x) x.alpha, CV_model, 'UniformOutput', true);
dsmax = cellfun(@(x) x.d, CV_model, 'UniformOutput', true);
Rmax = cellfun(@(x) x.R, CV_model, 'UniformOutput', true);
Nfea = cellfun(@(x) sum(abs(x.w)>0)-1, CV_model, 'UniformOutput', true);

figure(3)
    subplot(1,3,1)
    plot(alpha)
    xlim([0,101])
    subplot(1,3,2)
    plot(Rmax)
    xlim([0,101])
    subplot(1,3,3)
    plot(Nfea)
    xlim([0,101])
figure(4)
    subplot(1,3,1)
    surf(CVacc)
    ylim([0,101])
    subplot(1,3,2)
    surf(CVauc)
    ylim([0,101])
    subplot(1,3,3)
    surf(TCacc)
    ylim([0,101])
figure(3)
    subplot(1,2,1)
    surf(TCacc)
    ylim([0,101])
    subplot(1,2,2)
    surf(TPacc)
    ylim([0,101])




CV_acc = cell(10,1);
TC_acc = cell(10,1);
TP_acc = cell(10,1);
for i = 1:10
    load(sprintf('/storage/scratch/rr0311/0.hmp.output/HMPv13.c19.unifsR0.%02d.mat', i));
    CV_acc{i} = CVacc;
    TC_acc{i} = TCacc;
    TP_acc{i} = TPacc;
end
clearvars -except CV_acc TC_acc TP_acc;




figure(1)
hold on;
for k=1:10
    ca1 = cell2mat(cellfun(@(x) x(:,k)', CV_acc, 'UniformOutput', false))';
    ca2 = cell2mat(cellfun(@(x) x(:,k)', TC_acc, 'UniformOutput', false))';
    ca3 = cell2mat(cellfun(@(x) x(:,k)', TP_acc, 'UniformOutput', false))';
    ca4 = mean((ca2*250 + ca3*121)/371, 2);
    plot(ca4);
    ylim([0.5, 0.95]);
    pause(2)
end


figure(1)
    subplot(1,4,1)
    plot(ca1)
    xlim([0,110])
    subplot(1,4,2)
    plot(ca2)
    xlim([0,110])
    subplot(1,4,3)
    plot(ca3)
    xlim([0,110])
    subplot(1,4,4)
    plot(ca4)
    xlim([0,110])



%% Feature Selection Matrix %%
load('./0.hmp.input/HMPv13.phylo.mat')
load('./0.hmp.input/HMPv13.c19.mat')
load('./0.hmp.output/HMPv13.c19.unifsR2.mat')
subsample = sort([train.subsample;test.subsample]);
subtaxa = train.subtaxa(abs(CV_fs.w(2:end)) > 0);

fsmatrix = DataMatrix(HMPv13.edgebool.mat(:,subtaxa), HMPv13.edgebool.rownames, HMPv13.edgebool.colnames(subtaxa));
[~,J] = sort(sum(fsmatrix));
[~,I] = sort(fsmatrix.rownames);
fsmatrix = fsmatrix(I,J);
for k = size(fsmatrix,2):-1:1
    [~,I] = sort(double(fsmatrix(:,k)), 'descend');
    fsmatrix = fsmatrix(I,:);
end
fsmatrix = DataMatrix(full([CV_fs.w(abs(CV_fs.w)>0)';0,sum(fsmatrix);zeros(HMPv13.Ntaxa,1),double(fsmatrix)]), ['0Coefficient';'0Counts';fsmatrix.rownames], ['Intercept',fsmatrix.colnames]);
dmwrite(fsmatrix, sprintf('./0.hmp.table/HMPv13.c%02d.fs.table.%.4f.csv', [i,CV_fs.alpha]), 'Delimiter',',');

%% Distance Matrix %%
load('./0.hmp.input/HMPv13.phylo.mat')
load('./0.hmp.input/HMPv13.c19.mat')
load('./0.hmp.output/HMPv13.c19.unifsR.mat')
subsample = sort([train.subsample;test.subsample]);
subtaxa = train.subtaxa(abs(CV_fs.w(2:end)) > 0);

N = length(train.subsample);
edge_matrix.mat = HMPv13.edgematrix.mat(subsample, subtaxa);
edge_matrix.rownames = HMPv13.edgematrix.rownames(subsample);
edge_matrix.colnames = HMPv13.edgematrix.colnames(subtaxa);
edgedist = train.branch(abs(CV_fs.w(2:end)) > 0);
alpha = CV_fs.alpha;
R = CV_fs.d;
alpha0 = 1/(sqrt(N)*(1/alpha-1)+1);
R0 = alpha*R+sqrt(N)*(1-alpha)*R;

D0 = pdist(edge_matrix.mat, @(x,J) double((repmat(x,size(J,1),1)-J)).^2*edgedist);
dPCoA = DataMatrix(squareform(sqrt(D0)), edge_matrix.rownames, edge_matrix.rownames);
dmwrite(dPCoA, './0.hmp.table/HMPv13.c19.fs.dist.dPCoA.csv', 'Delimiter',',');
D1 = pdist(edge_matrix.mat, @(x,J) abs(double(repmat(x,size(J,1),1)-J))*edgedist);
unifrac = DataMatrix(squareform(D1), edge_matrix.rownames, edge_matrix.rownames);
dmwrite(unifrac, './0.hmp.table/HMPv13.c19.fs.dist.w.non_nor.csv', 'Delimiter',',');
Dalpha = pdist(edge_matrix.mat*100, @(x,J) arrayfun(@(y) unifsR_dis(x,J(y,:),edgedist,edgedist,alpha,N), 1:size(J,1), 'UniformOutput', true));
mydist = DataMatrix(squareform(Dalpha), edge_matrix.rownames, edge_matrix.rownames);
dmwrite(mydist, './0.hmp.table/HMPv13.c19.fs.dist.alpha.csv', 'Delimiter',',');
Dfs = alpha0/R0*D1 + (1-alpha0)/R0^2*D0;
mydist = DataMatrix(squareform(Dfs), edge_matrix.rownames, edge_matrix.rownames);
dmwrite(mydist, './0.hmp.table/HMPv13.c19.fs.dist.comb.csv', 'Delimiter',',');




HMPv13.tax_table.Tax_ID = HMPv13.tax_table.Tax_ID';

taxalist = struct2cell(HMPv13.tax_table)';
taxalist = cat(2,taxalist{:});
brlist = DataMatrix(HMPv13.edgebool.mat, HMPv13.edgebool.rownames, HMPv13.edgebool.colnames);
brlist = brlist(HMPv13.tax_table.Tax_ID, :);
cutoff = 50;

parfor k = 1:size(brlist,2)
    rm_node(k) = 0;
    index = brlist(:,k)==1;
    infolist = taxalist(index,2:end-1);
    infolist = [infolist, arrayfun(@(x) strjoin(infolist(x,:),''), 1:sum(index), 'UniformOutput',false)'];
    if sum(index)<=cutoff
        [fre,~,~,lab] = crosstab(infolist(:,end));
        if any(fre/sum(fre)>=0.8)
            infolist = infolist(strcmp(infolist(:,end),lab(fre/sum(fre)>=0.8)),:);
            rm_node(k) = 1;
            newlist{k} = [brlist.colnames(k),infolist(1,1:5)];
        end
    end
end
tax_table = cat(1,newlist{:});
fid = fopen('file.csv','wt');
for i=1:size(tax_table,1)
    fprintf(fid, '%s,%s,%s,%s,%s,%s\n', tax_table{i,:});
end
fclose(fid);

rm_node(train.subtaxa) = 0;
newtree = prune(HMPv13.phy_tree, logical(rm_node), 'Mode','Exclusive');
phytreewrite('./0.hmp.table/HMPv13.c19.ref.tre', newtree);




taxalist = struct2cell(phylo.tax_table)';
taxalist = cat(2,taxalist{:});
brlist = DataMatrix(phylo.edgebool.mat, phylo.edgebool.rownames, phylo.edgebool.colnames);
brlist = brlist(phylo.tax_table.Tax_ID, :);
cutoff = 5;

parfor k = 1:size(brlist,2)
    index = brlist(:,k)==1;
    otulist = phylo.tax_table.Order(index);
    rm_node(k) = 0;
    if sum(index)<=cutoff && all(strcmp(otulist(1),otulist))
        otuinfo = taxalist(index,:);
        rm_node(k) = 1;
        newlist{k} = [brlist.colnames(k),otuinfo(1,2:7)];
    end
end
newlist = cat(1,newlist{:});
newtree = prune(phylo.phy_tree, logical(rm_node), 'Mode','Exclusive');











[1:numLeaves] for the leaves and as [numLeaves+1:numLeaves+numBranches]


addpath(genpath(pwd));
import bioma.data.*
clear;
load('./0.matlab.out/oral_new.phylo.mat');
load('./0.matlab.out/oral_new.unifsR.mat')

dmwrite(DataMatrix(full(HMPv13.edgematrix.mat*100), phylo.edgematrix.rownames, phylo.edgematrix.colnames), './0.matlab.out/oral_new.edgematrix.csv', 'Delimiter', ',');
phytreewrite('./0.matlab.out/oral_new.ref.tre', phylo.phy_tree);
csvwrite('./0.matlab.out/oral_new.subtaxa.csv', train.subtaxa');




dmwrite(DataMatrix(coef1, phylo_train.edgebool.colnames, parlist.alpha), './0.matlab.out/HMP_100_Att_Sup.unifs_cv1.weight.csv', 'Delimiter', ',');

ta = cellfun(@(x) confusionmat(test.bodysite,nominal(x)), CV_pred, 'UniformOutput', false);
tb = cellfun(@(x) (x(1)+x(4))/sum(x(:)), ta, 'UniformOutput', true);


tic;
[CV_model, CV_score, CV_opt, time0] = unifsR_cv(phylo_train.edgematrix.mat*100, phylo_train.sam_table.HMPbodysubsite, get(phylo_train.phy_tree, 'Distance'), []);
time = toc;

unifs_cv1 = unifsnew1_crossval(phylo_train.edgematrix.mat, phylo_train.sam_table.HMPbodysubsite, get(phylo_train.phy_tree, 'Distance'), 10, parlist);
time1 = toc;
disp(time1);
coef1 = cell2mat(cellfun(@(x) x.cvmodel.w, unifs_cv1, 'UniformOutput', false));
[pred1, score1] = cellfun(@(x) unifs_pred(x.cvmodel, phylo_test), unifs_cv1, 'UniformOutput', false);
confmat1 = cellfun(@(x) confusionmat(cellstr(x), phylo_test.sam_table.HMPbodysubsite), pred1, 'UniformOutput', false);
disp(sum(cat(3,confmat1{:}),3));
[X1, Y1, T1, AUC1] = cellfun(@(x) perfcurve(phylo_test.sam_table.HMPbodysubsite, x, 'Supragingival_plaque'), score1, 'UniformOutput', false);
score = DataMatrix([cell2mat(cellfun(@(x) x.score, unifs_cv1, 'UniformOutput', false));cell2mat(score1)],[phylo_train.sam_table.Sam_ID;phylo_test.sam_table.Sam_ID], parlist.alpha);
save('./0.matlab.out/HMP_100_Att_Sup.unifs_cv1.mat', 'unifs_cv1', 'time1', 'pred1', 'score1', 'confmat1', 'X1', 'Y1', 'T1','AUC1');
dmwrite(DataMatrix(coef1, phylo_train.edgebool.colnames, parlist.alpha), './0.matlab.out/HMP_100_Att_Sup.unifs_cv1.weight.csv', 'Delimiter', ',');
dmwrite(score, './0.matlab.out/HMP_100_Att_Sup.unifs_cv1.score.csv', 'Delimiter', ',');




coef = zeros(Nexd+1, kalpha);
coef(1+[0,train.subtaxa], :) = cell2mat(cellfun(@(x) x.model.w, CV_alpha', 'UniformOutput', false));
coef = DataMatrix(coef, ['Intercept';HMPv13.edgematrix.colnames], cellfun(@(x) x.model.alpha, CV_alpha, 'UniformOutput', true));
dmwrite(coef, './0.hmp.output/HMPv13.c19.CV_alpha.weight.csv', 'Delimiter', ',');

coef = zeros(Nexd+1, kD);
coef(1+[0,train.subtaxa], :) = cell2mat(cellfun(@(x) x.model.w, CV_D, 'UniformOutput', false));
coef = DataMatrix(coef, ['Intercept';HMPv13.edgematrix.colnames], cellfun(@(x) x.model.D, CV_D, 'UniformOutput', true));
dmwrite(coef, './0.hmp.output/HMPv13.c43.CV_D.weight.csv', 'Delimiter', ',');

ft = zeros(Nexd+1, kalpha);
ft(1+[0,train.subtaxa], :) = cell2mat(cellfun(@(x) x.w, CV_model(:,10)', 'UniformOutput', false));
ft = ft(2:end,end); 

otulist = HMPv13.edgebool.mat(:,ft~=0);
rowname = HMPv13.edgebool.rownames;
colname = HMPv13.edgebool.colnames(ft~=0);
otulist = DataMatrix(full(otulist), rowname, colname);
dmwrite(otulist, './0.hmp.output/HMPv13.c19.otulist.csv', 'Delimiter', ',');


score = DataMatrix([cell2mat(cellfun(@(x) x.score, CV_alpha, 'UniformOutput', false));cell2mat(score1)],[phylo_train.sam_table.Sam_ID;phylo_test.sam_table.Sam_ID], parlist.alpha);

save('./0.matlab.out/HMP_100_Att_Sup.unifs_cv1.mat', 'unifs_cv1', 'time1', 'pred1', 'score1', 'confmat1', 'X1', 'Y1', 'T1','AUC1');
dmwrite(DataMatrix(coef1, phylo_train.edgebool.colnames, parlist.alpha), './0.matlab.out/HMP_100_Att_Sup.unifs_cv1.weight.csv', 'Delimiter', ',');
dmwrite(score, './0.matlab.out/HMP_100_Att_Sup.unifs_cv1.score.csv', 'Delimiter', ',');


alpha0 = 1/(sqrt(200)*(1/alpha-1)+1);

CV_alpha{19}.model.alpha

mat = [train.edgematrix;test.edgematrix];
index = CV_alpha{19}.model.w~=0;
mat = mat(:,index(2:end));
edgedist = train.branch(index(2:end));


[~,~,edgematrix,distance] = fastUnifrac(HMPv13.otu_table, HMPv13.phy_tree, [], {'edge','dPCoA','non_w','w.non_nor','w.nor'});
[~,~,edgematrix,distance] = fastUnifrac(phylo.otu_table, phylo.phy_tree, [], {'edge','dPCoA','non_w','w.non_nor','w.nor'});


edge_matrix.mat = HMPv13.edgematrix.mat;
edge_matrix.rownames = HMPv13.edgematrix.rownames;
edgedist = get(HMPv13.phy_tree, 'Distance');

D = pdist(edge_matrix.mat, @(x,J) sqrt(double((repmat(x,size(J,1),1)-J)).^2*edgedist));
dist_list.dPCoA = DataMatrix(squareform(D), edge_matrix.rownames, edge_matrix.rownames);
dmwrite(dist_list.dPCoA, './0.hmp.output/HMPv13.dist.dPCoA.csv', 'Delimiter',',');

D = pdist(edge_matrix.mat, @(x,J) (abs(double(repmat(x>0,size(J,1),1)-(J>0)))*edgedist)./(abs(double(repmat(x>0,size(J,1),1)|(J>0)))*edgedist));
dist_list.non_w = DataMatrix(squareform(D), edge_matrix.rownames, edge_matrix.rownames);
dmwrite(dist_list.non_w, './0.hmp.output/HMPv13.dist.non_w.csv', 'Delimiter',',');

D = pdist(edge_matrix.mat, @(x,J) abs(double(repmat(x,size(J,1),1)-J))*edgedist);
dist_list.w.non_nor = DataMatrix(squareform(D), edge_matrix.rownames, edge_matrix.rownames);
dmwrite(dist_list.w.non_nor, './0.hmp.output/HMPv13.dist.w.non_nor.csv', 'Delimiter',',');

D = pdist(edge_matrix.mat, @(x,J) (abs(double(repmat(x,size(J,1),1)-J))*edgedist)./(abs(double(repmat(x,size(J,1),1)+J))*edgedist));
dist_list.w.nor = DataMatrix(squareform(D), edge_matrix.rownames, edge_matrix.rownames);
dmwrite(dist_list.w.nor, './0.hmp.output/HMPv13.dist.w.nor.csv', 'Delimiter',',');

parfor k = 1:10
    [~,~,b] = intersect(unifs_cv{6,k}{4}(:,end), unifs_cv{6,end}{4}(:,end));
    cpp{k} = b;
end









x = [100./linspace(90,20,30)];
plot(x, arrayfun(@(x) CV_model{101,x}.delta*CV_model{101,x}.D(2), 1:30, 'UniformOutput', true), 'Color', 'red', 'LineStyle', ':');
hold
plot(x, arrayfun(@(x) CV_model{1,x}.delta*CV_model{1,x}.D(2), 1:30, 'UniformOutput', true), 'Color', 'blue', 'LineStyle', ':');


[CV_fs.index/101, CVacc(CV_fs.index), max(CVacc(:)), TPacc(CV_fs.index), max(TPacc(:))]



plot(x, arrayfun(@(x) CV_model{1,x}.delta*CV_model{1,x}.D(2), 1:30, 'UniformOutput', true), 'Color', 'blue', 'LineStyle', ':');



plotC = cell(size(unifs_cv,1),1);
plotD = cell(size(unifs_cv,1),1);
plotR = cell(size(unifs_cv,1),1);
plotalpha = cell(size(unifs_cv,1),1);
plotlambda = cell(size(unifs_cv,1),1);
plotdelta = cell(size(unifs_cv,1),1);
for k = 1:size(unifs_cv,1)
    ck = cellfun(@(x) x{6}, unifs_cv(k,:), 'UniformOutput', false);
    ck = cat(2,ck{:});
    plotD{k} = cellfun(@(x) x.D(2), ck, 'UniformOutput', true);
    plotdelta{k} = cellfun(@(x) x.delta, ck, 'UniformOutput', true);
    plotR{k} = cellfun(@(x) x.R, ck, 'UniformOutput', true);
    plotalpha{k} = cellfun(@(x) x.alpha, ck, 'UniformOutput', true);
    plotlambda{k} = cellfun(@(x) x.R*x.alpha, ck, 'UniformOutput', true);
    plotC{k} = cellfun(@(x) x.D(2)*x.delta, ck, 'UniformOutput', true);
end

[CV_fs.index/101, CVacc(CV_fs.index), max(CVacc(:)), TPacc(CV_fs.index), max(TPacc(:))]

mean(plotC{2}(end,1:10))


k=9


x = cellfun(@(x) x(1,k), plotD, 'UniformOutput', true);
y1 = cellfun(@(x) x(1,k), plotC, 'UniformOutput', true);
y2 = cellfun(@(x) x(end,k), plotC, 'UniformOutput', true);
plot(x, y1, 'Color','blue', 'LineStyle', ':')
plot(x,y2, 'Color','red', 'LineStyle', ':')



