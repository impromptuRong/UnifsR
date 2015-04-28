%% Simulation Cross Validation Test: 20 samples each, 100 OTUs, 10 G-Featres, 5 G-Noises  %%
clear;
% Read In simulated Probability %
% view(reroot(phytreeread('sim.100.tre')))
prob_table = DataMatrix('File', 'sim.prob.100.csv', 'Delimiter', ',');
rownames = regexprep(strrep(prob_table.RowNames, '"', ''), '[^a-zA-Z0-9]', '_');
colnames = regexprep(strrep(prob_table.ColNames, '"', ''), '[^a-zA-Z0-9]', '_');
shuff_list = arrayfun(@(i) tabulate(randsample(1:size(prob_table,2), 1000, true, double(prob_table(i,:)))), 1:size(prob_table,1), 'UniformOutput', false);
rand_table = cell2mat(cellfun(@(k) [k(:,2)', zeros(1,size(prob_table,2)-size(k,1))], shuff_list', 'UniformOutput', false));

% No noise and outlier %
otu_table = DataMatrix(rand_table, 'RowNames', rownames, 'ColNames', colnames);
dmwrite(otu_table, './0.matlab.out/sim_100.summary.csv', 'Delimiter', ',');
phylo = phyloseq('otu_table', 'sim_100.summary.csv', 'phy_tree', 'sim.100.tre', 'sam_table', 'sim.info.100.csv');
save('./0.matlab.out/sim_100.phylo.mat', 'phylo');
phytreewrite('./0.matlab.out/sim_100.ref.tre', phylo.phy_tree);
dmwrite(DataMatrix(full(phylo.edgebool.mat), phylo.edgebool.rownames, phylo.edgebool.colnames), './0.matlab.out/sim_100.ref.tre.edge.csv', 'Delimiter', ',');

clear;
parlist.alpha = 0:0.1:1;
parlist.D = 1./[1:5;1:5];
load('./0.matlab.out/sim_100.phylo.mat');

tic;
unifs_cv1 = unifs_crossval1(phylo.edgematrix.mat, phylo.sam_table.Group, get(phylo.phy_tree, 'Distance'), 10, parlist);
time1 = toc;
disp(time1);
save('./0.matlab.out/sim_100.unifs_cv1.mat', 'unifs_cv1', 'time1');
coef1 = cell2mat(cellfun(@(x) x.cvmodel.w, unifs_cv1, 'UniformOutput', false));
dmwrite(DataMatrix(coef1, phylo.edgebool.colnames, parlist.alpha), './0.matlab.out/sim_100.unifs_cv1.weight.csv', 'Delimiter', ',');

tic;
unifs_cv2 = unifs_crossval2(phylo.edgematrix.mat, phylo.sam_table.Group, get(phylo.phy_tree, 'Distance'), 10, parlist);
time2 = toc;
disp(time2);
save('./0.matlab.out/sim_100.unifs_cv2.mat', 'unifs_cv2', 'time2');
coef2 = cell2mat(cellfun(@(x) x.cvmodel.w, unifs_cv2, 'UniformOutput', false));
dmwrite(DataMatrix(coef2, phylo.edgebool.colnames, parlist.alpha), './0.matlab.out/sim_100.unifs_cv2.weight.csv', 'Delimiter', ',');

% Introduce noise and outlier %
clear;
rand_table = DataMatrix('File', 'sim_100.summary.csv', 'Delimiter', ',');
rownames = regexprep(strrep(rand_table.RowNames, '"', ''), '[^a-zA-Z0-9]', '_');
colnames = regexprep(strrep(rand_table.ColNames, '"', ''), '[^a-zA-Z0-9]', '_');
outlier = ceil((1000.^(sprand(size(rand_table,1),size(rand_table,2),0.1))-1)/10);
noise = ceil(sprand(size(rand_table,1),size(rand_table,2),1)*20);
noise = full((noise-10).*(noise>0));
otu_table = double(rand_table) + outlier + noise;
otu_table = otu_table.*(otu_table>=0);
%     noise = sprand(phylo.Nsample, phylo.Ntaxa, 1);
%     outlier = sqrt(sprand(phylo.Nsample,phylo.Ntaxa,0.01)+1);
%     adj = ceil(log(sprand(phylo.Nsample,phylo.Ntaxa,0.01)+1)*100) + (ceil(noise*80)-40).*(noise>0);
%     spy(otu_table);
%     nnz(otu_table)/prod(size(otu_table))
% Construct Phyloseq Class %
otu_table = DataMatrix(otu_table, 'RowNames', rownames, 'ColNames', colnames);
dmwrite(otu_table, './0.matlab.out/sim_100_sp.summary.csv', 'Delimiter', ',');
phylo = phyloseq('otu_table', 'sim_100_sp.summary.csv', 'phy_tree', 'sim.100.tre', 'sam_table', 'sim.info.100.csv');
save('./0.matlab.out/sim_100_sp.phylo.mat', 'phylo');
phytreewrite('./0.matlab.out/sim_100.ref.tre', phylo.phy_tree);
dmwrite(DataMatrix(full(phylo.edgebool.mat), phylo.edgebool.rownames, phylo.edgebool.colnames), './0.matlab.out/sim_100.ref.tre.edge.csv', 'Delimiter', ',');

clear;
parlist.alpha = 0:0.1:1;
parlist.D = 1./[1:5;1:5];
load('./0.matlab.out/sim_100_sp.phylo.mat');

tic;
unifs_cv1 = unifs_crossval1(phylo.edgematrix.mat, phylo.sam_table.Group, get(phylo.phy_tree, 'Distance'), 10, parlist);
time1 = toc;
disp(time1);
save('./0.matlab.out/sim_100_sp.unifs_cv1.mat', 'unifs_cv1', 'time1');
coef1 = cell2mat(cellfun(@(x) x.cvmodel.w, unifs_cv1, 'UniformOutput', false));
dmwrite(DataMatrix(coef1, phylo.edgebool.colnames, parlist.alpha), './0.matlab.out/sim_100_sp.unifs_cv1.weight.csv', 'Delimiter', ',');

tic;
unifs_cv2 = unifs_crossval2(phylo.edgematrix.mat, phylo.sam_table.Group, get(phylo.phy_tree, 'Distance'), 10, parlist);
time2 = toc;
disp(time2);
save('./0.matlab.out/sim_100_sp.unifs_cv2.mat', 'unifs_cv2', 'time2');
coef2 = cell2mat(cellfun(@(x) x.cvmodel.w, unifs_cv2, 'UniformOutput', false));
dmwrite(DataMatrix(coef2, phylo.edgebool.colnames, parlist.alpha), './0.matlab.out/sim_100_sp.unifs_cv2.weight.csv', 'Delimiter', ',');

