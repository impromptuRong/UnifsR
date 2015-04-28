function phylo = subset_phylo(phylo, subsample, subtaxa)
%% Initialize Input Variables %%
import bioma.data.*
if(isempty(subsample))
    subsample = 1:phylo.Nsample;
end
if(isempty(subtaxa))
    subtaxa = 1:phylo.Ntaxa;
end

phylo.sam_table = structfun(@(x) x(subsample), phylo.sam_table, 'UniformOutput', false);
phylo.tax_table = structfun(@(x) x(subtaxa), phylo.tax_table, 'UniformOutput', false);
phylo.otu_table = phylo.otu_table(subsample, subtaxa);
tmp = 1:phylo.Ntaxa;
nodename = get(phylo.phy_tree, 'NodeNames');
phylo.phy_tree = prune(phylo.phy_tree, setdiff(tmp, tmp(subtaxa)));
subnode = ismember(nodename, get(phylo.phy_tree, 'NodeNames'));
phylo.edgematrix.mat = phylo.edgematrix.mat(subsample, subnode);
phylo.edgematrix.rownames = phylo.edgematrix.rownames(subsample);
phylo.edgematrix.colnames = phylo.edgematrix.colnames(subnode);
phylo.edgebool.mat = phylo.edgebool.mat(subtaxa, subnode);
phylo.edgebool.rownames = phylo.edgebool.rownames(subtaxa);
phylo.edgebool.colnames = phylo.edgebool.colnames(subnode);
if(isfield(phylo, 'edgelist'))
    phylo.edgelist = phylo.edgelist(subnode);
end
[phylo.Nsample, phylo.Ntaxa] = size(phylo.otu_table);
phylo.Nexd = size(phylo.edgematrix.mat,2);
end
