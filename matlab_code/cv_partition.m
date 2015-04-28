function [ps] = cv_partition(phylo, class, label, parlist)
%% Initialize Input Variables %%
if(~isfield(parlist,'trim') || isempty(parlist.trim))
    parlist.trim = true;
end
if(islogical(parlist.trim) && parlist.trim)
    parlist.trim = [0; 0; 0];
end

%% Extract Samples with Label defined in Class %%
seqdep = phylo.seqdep;
high_abu = seqdep > parlist.trim(3);
group = phylo.sam_table.(class);
if(isempty(label))
    label = unique(group(high_abu));
end
nlevel = length(label);
subindex = sum(cell2mat(arrayfun(@(x) strcmp(group, label{x})*x, 1:nlevel, 'UniformOutput', false)),2);
subtrain = high_abu.*subindex;
taxa_tab = double(phylo.otu_table./repmat(seqdep,1,phylo.Ntaxa));
if(~islogical(parlist.trim))
    subtaxa = find(mean(taxa_tab(subtrain~=0,:),1) > parlist.trim(1));
else
    subtaxa = 1:phylo.Ntaxa;
end
if(isfield(phylo, 'edgematrix') && ~isempty(phylo.edgematrix.mat))
    node_tab = phylo.edgematrix.mat;
    if(~islogical(parlist.trim))
        subnode = find(mean(node_tab(subtrain~=0,:),1) > parlist.trim(2));
    else
        subnode = 1:phylo.Nexd;
    end
end

%% Parameters for Splitting and CV %%
nclass = crosstab(subindex(subindex~=0));
if(~isfield(parlist,'size') || isempty(parlist.size))
    parlist.size = 2/3;
end
if(length(parlist.size)==1)
    parlist.size = parlist.size*ones(nlevel,1);
end
train_size = parlist.size;
train_size(train_size<=1) = train_size(train_size<=1).*nclass(train_size<=1);
train_size = ceil(train_size);
if(~isfield(parlist,'kfold') || isempty(parlist.kfold))
    parlist.kfold = 10;
end

%% Split Training and Testing, Stratified CV id %%
train_label = sort(cell2mat(arrayfun(@(x) randsample(find(subtrain==x), train_size(x))', 1:nlevel, 'UniformOutput', false)))';
test_label = sort(setdiff(find(subindex~=0), train_label));
train_group = nominal(group(train_label));
test_group = nominal(group(test_label));
kfold = parlist.kfold;
cp = cvpartition(train_group, 'k', kfold);
foldid = sum(cell2mat(arrayfun(@(x) cp.test(x)*x, 1:cp.NumTestSets, 'UniformOutput',false)),2);

%% Output Partition Structure %%
ps = struct();
ps.train_label = train_label;
ps.train_group = train_group;
ps.test_label = test_label;
ps.test_group = test_group;
ps.foldid = foldid;
if(isfield(phylo, 'edgematrix') && ~isempty(phylo.edgematrix.mat))
    ps.subnode = subnode;
end
ps.subtaxa = subtaxa;
ps.size = [length(train_label)/(length(train_label)+length(test_label)); kfold];
ps.trim = parlist.trim;
end
