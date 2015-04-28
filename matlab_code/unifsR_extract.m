addpath(genpath(pwd));
import bioma.data.*
load('./0.hmp.input/HMPv13.phylo.mat')
load('./0.hmp.input/HMPv13.c19.00.mat')
load('./0.hmp.output/HMPv13.c19.unifsR1.00.mat')
import bioma.data.*
% subsample = sort([train.subsample;test.subsample]);
% subtaxa = train.subtaxa(abs(CV_fs.w(2:end)) > 0);

phylo0 = HMPv13;
model0 = CV_fs;
model0.subtaxa = train.subtaxa;
model0.subsample = sort([train.subsample;test.subsample]);
slot0 = {'HMPbodysubsite', 'genus'};

clearvars -except model0 phylo0 slot0
phylo = phylo0;
model = model0;
slot = slot0;

fsmatrix = DataMatrix(edges.mat(:,subtaxa), HMPv13.edgebool.rownames, HMPv13.edgebool.colnames(subtaxa));
fsmatrix = edges;
tax_table = tax_table0;

function [] = unifsR_extract(model, phylo, slot)
%% Feature Selection Matrix %%
if ~isfield(model, 'subtaxa')
    model.subtaxa = 1:length(model.w)-1;
    endq
model.subtaxa = model.subtaxa(abs(model.w(2:end)) > 0);
fsmatrix = DataMatrix(double(phylo.edgebool.mat(:,model.subtaxa)), phylo.edgebool.rownames, phylo.edgebool.colnames(model.subtaxa));
otu_table = phylo.otu_table./repmat(sum(phylo.otu_table,2), 1, phylo.Ntaxa)*100;
otu_table = otu_table(model.subsample,:);
group = phylo.sam_table.(slot{1})(model.subsample);
g1 = strcmp(group, model.label{1});
g2 = strcmp(group, model.label{2});
trownames = {[model.label{1},'.mean'], [model.label{1},'.pc'], [model.label{1},'.se'], [model.label{2},'.mean'], [model.label{2},'.pc'], [model.label{2},'.se'], 'Ttest', 'Wilcox'};
tax_table = phylo.tax_table;
% clearvars -except model fsmatrix otu_table group tax_table slot;

%% Organize Features with Enterotype %%
[Notu, Nfea] = size(fsmatrix);
[~, J] = sort(sum(fsmatrix));
fsmatrix = -fsmatrix(:,J);
entype = cell(Nfea, 3);
statest = cell(Nfea, 3);
fsummary = zeros(Nfea, 12);

%% Statistics for Enterotypes and Features %%
for k = 1:Nfea
    %% Enterotype and Features List %%
    subotu_ID = fsmatrix.rownames(fsmatrix(:,k)==-1);
    fsmatrix(subotu_ID, k:Nfea) = -fsmatrix(subotu_ID, k:Nfea)*k;
    subfea_ID = unique(full(double(fsmatrix(fsmatrix(:,k)>0, k))));
    entype(k, 1:2) = {subotu_ID, subfea_ID};
    
    % Enterotype: Single OTU under %
    sub_table = double(otu_table(:, subotu_ID));
    sub_total = [sum(sum(sub_table(g1,:))), sum(sum(sub_table(g2,:)))];
    parfor m = 1:length(subotu_ID)
        [~, p1] = ttest2(sub_table(g1,m), sub_table(g2,m), 'Vartype', 'unequal');
        [p2, ~] = ranksum(sub_table(g1,m), sub_table(g2,m));
        output1(m,:) = [mean(sub_table(g1,m)), sum(sub_table(g1,m))/sub_total(1), std(sub_table(g1,m))/sqrt(sum(g1)), ...
            mean(sub_table(g2,m)), sum(sub_table(g2,m))/sub_total(2), std(sub_table(g2,m))/sqrt(sum(g2)), p1, p2];
    end
    output1 = DataMatrix(output1, subotu_ID, trownames);
    
    % Enterotype: Combined Level %
    [~, I] = ismember(subotu_ID, tax_table.Tax_ID);
    sub_taxa = tax_table.(slot{2})(I);
    sub_table = cellfun(@(x) sum(sub_table(:,strcmp(sub_taxa, x)),2), unique(sub_taxa), 'UniformOutput', false);
    sub_total = [sum(cellfun(@(x) sum(x(g1)), sub_table, 'UniformOutput',true)), sum(cellfun(@(x) sum(x(g2)), sub_table, 'UniformOutput',true))];
    parfor m = 1:length(sub_table)
        [~, p1] = ttest2(sub_table{m}(g1), sub_table{m}(g2), 'Vartype', 'unequal');
        [p2, ~] = ranksum(sub_table{m}(g1), sub_table{m}(g2));
        output2(m,:) = [mean(sub_table{m}(g1)), sum(sub_table{m}(g1))/sub_total(1), std(sub_table{m}(g1))/sqrt(sum(g1)), ...
            mean(sub_table{m}(g2)), sum(sub_table{m}(g2))/sub_total(2), std(sub_table{m}(g2))/sqrt(sum(g2)), p1, p2];
    end
    output2 = DataMatrix(output2, unique(sub_taxa), trownames);
    
    % Enterotype: Whole Enterotype %
    sub_table = sum(cell2mat(sub_table'), 2);
    [~, p1] = ttest2(sub_table(g1), sub_table(g2), 'Vartype', 'unequal');
    [p2, ~] = ranksum(sub_table(g1), sub_table(g2));
    output3 = [mean(sub_table(g1)), std(sub_table(g1))/sqrt(sum(g1)), mean(sub_table(g2)), std(sub_table(g2))/sqrt(sum(g2)), p1, p2];
    
    % Feature list: Combined Level %
    subfea_ID = entype(subfea_ID, 1);
    subfea_ID = cat(1, subfea_ID{:});
    [~, I] = ismember(subfea_ID, tax_table.Tax_ID);
    sub_taxa = tax_table.(slot{2})(I);
    sub_table = double(otu_table(:, subfea_ID));
    sub_table = cellfun(@(x) sum(sub_table(:,strcmp(sub_taxa, x)),2), unique(sub_taxa), 'UniformOutput', false);
    sub_total = [sum(cellfun(@(x) sum(x(g1)), sub_table, 'UniformOutput',true)), sum(cellfun(@(x) sum(x(g2)), sub_table, 'UniformOutput',true))];
    parfor m = 1:length(sub_table)
        [~, p1] = ttest2(sub_table{m}(g1), sub_table{m}(g2), 'Vartype', 'unequal');
        [p2, ~] = ranksum(sub_table{m}(g1), sub_table{m}(g2));
        output4(m,:) = [mean(sub_table{m}(g1)), sum(sub_table{m}(g1))/sub_total(1), std(sub_table{m}(g1))/sqrt(sum(g1)), ...
            mean(sub_table{m}(g2)), sum(sub_table{m}(g2))/sub_total(2), std(sub_table{m}(g2))/sqrt(sum(g2)), p1, p2];
    end
    output4 = DataMatrix(output3, unique(sub_taxa), trownames);
    
    % Feature list: Whole Enterotype %
    sub_table = sum(cell2mat(sub_table'), 2);
    [~, p1] = ttest2(sub_table(g1), sub_table(g2), 'Vartype', 'unequal');
    [p2, ~] = ranksum(sub_table(g1), sub_table(g2));
    output5 = [mean(sub_table(g1)), std(sub_table(g1))/sqrt(sum(g1)), mean(sub_table(g2)), std(sub_table(g2))/sqrt(sum(g2)), p1, p2];
    
    % Store Result and eliminate values %
    statest(k, :) = {output1, output2, output3, output4, output5};
    clearvars output1 output2 output3 output4 output5;
end
fsummary = DataMatrix(fsummary, fsmatrix.colnames, [trownames([1,3,4,6,7,8]),trownames([1,3,4,6,7,8])]);

entype(:,3) = {[]};
subotu_ID = 0;
for k = 1:Nfea
    if(isempty(entype{k,3}))
        subotu_ID = subotu_ID+1;
        entype(k, 3) = {subotu_ID};
        ref_abu = statest{k,2}(:,[2,5]);
        
        cand_list = find(cellfun(@(x) ismember(k,x), entype(k+1:end,2), 'UniformOutput', true)) + k;
        common = k:Nfea;
        for m = 1:length(cand_list)
            common = intersect(common, entype{cand_list(m),2});
        end
        common = setdiff(common, find(cellfun(@(x) ~isempty(x), entype(:,3))));
        
        for m = 1:length(common)
            com_abu = statest{common(m),2}(:,[2,5]);
            tar_abu = zeros(ref_abu.NRows, 2);
            for r = 1:ref_abu.NRows
                if ismember(ref_abu.rownames(r), com_abu.rownames)
                    tar_abu(r, :) = double(com_abu(ref_abu.rownames(r), :));
                end
            end
            if(all(sum(abs(ref_abu-tar_abu))<0.15))
                entype(common(m), 3) = {subotu_ID};
            end
        end
    end
end

entype2 = cell(entype_ID, 5);
for k = 1:entype_ID
    subfea_ID = find(cell2mat(entype(:,3))==k);
    subotu_ID = cat(1,entype{subfea_ID,1});
    
    % Enterotype: Single OTU %
    sub_table = double(otu_table(:, subotu_ID));
    sub_total = [sum(sum(sub_table(g1,:))), sum(sum(sub_table(g2,:)))];
    parfor m = 1:length(subotu_ID)
        [~, p1] = ttest2(sub_table(g1,m), sub_table(g2,m), 'Vartype', 'unequal');
        [p2, ~] = ranksum(sub_table(g1,m), sub_table(g2,m));
        output1(m,:) = [mean(sub_table(g1,m)), sum(sub_table(g1,m))/sub_total(1), std(sub_table(g1,m))/sqrt(sum(g1)), ...
            mean(sub_table(g2,m)), sum(sub_table(g2,m))/sub_total(2), std(sub_table(g2,m))/sqrt(sum(g2)), p1, p2];
    end
    output1 = DataMatrix(output1, subotu_ID, trownames);
    
    % Enterotype: Combined Level %
    [~, I] = ismember(subotu_ID, tax_table.Tax_ID);
    sub_taxa = tax_table.(slot{2})(I);
    sub_table = cellfun(@(x) sum(sub_table(:,strcmp(sub_taxa, x)),2), unique(sub_taxa), 'UniformOutput', false);
    sub_total = [sum(cellfun(@(x) sum(x(g1)), sub_table, 'UniformOutput',true)), sum(cellfun(@(x) sum(x(g2)), sub_table, 'UniformOutput',true))];
    parfor m = 1:length(sub_table)
        [~, p1] = ttest2(sub_table{m}(g1), sub_table{m}(g2), 'Vartype', 'unequal');
        [p2, ~] = ranksum(sub_table{m}(g1), sub_table{m}(g2));
        output2(m,:) = [mean(sub_table{m}(g1)), sum(sub_table{m}(g1))/sub_total(1), std(sub_table{m}(g1))/sqrt(sum(g1)), ...
            mean(sub_table{m}(g2)), sum(sub_table{m}(g2))/sub_total(2), std(sub_table{m}(g2))/sqrt(sum(g2)), p1, p2];
    end
    output2 = DataMatrix(output2, unique(sub_taxa), trownames);
    
    % Enterotype: Whole Enterotype %
    sub_table = sum(cell2mat(sub_table'), 2);
    [~, p1] = ttest2(sub_table(g1), sub_table(g2), 'Vartype', 'unequal');
    [p2, ~] = ranksum(sub_table(g1), sub_table(g2));
    output3 = [mean(sub_table(g1)), std(sub_table(g1))/sqrt(sum(g1)), mean(sub_table(g2)), std(sub_table(g2))/sqrt(sum(g2)), p1, p2];
    
    entype2(k,:) = {subfea_ID, subotu_ID, output3, output1, output2};
    clearvars output1 output2 output3;
end





%% Combine Similar Enterotype %%

%     field = fields(tax_table);
%     for taxa = sub_taxa
%         tax_table.(field{k}) = tax_table.(field{k})(I);
%     end

% for k = Nfea:-1:1
%     [~, I] = sort(double(fsmatrix(:,k)), 'descend');
%     fsmatrix = fsmatrix(I,:);
% end

%% OTU Matrix %%
group = phylo.sam_table.(slot{1})(model.subsample);
[~, I] = ismember(fsmatrix.rownames, tax_table.Tax_ID);
otu_table = phylo.otu_table(model.subsample, I);
field = fields(tax_table);
for k = 1:numel(field)
    tax_table.(field{k}) = tax_table.(field{k})(I);
end

fsmatrix = DataMatrix(full([CV_fs.w(abs(CV_fs.w)>0)';0,sum(fsmatrix);zeros(HMPv13.Ntaxa,1),double(fsmatrix)]), ['0Coefficient';'0Counts';fsmatrix.rownames], ['Intercept',fsmatrix.colnames]);
dmwrite(fsmatrix, './0.hmp.table/try.csv', 'Delimiter',',');

DataMatrix(fsmatrix.mat, fsmatrix.rownames, fsmatrix.colnames);







