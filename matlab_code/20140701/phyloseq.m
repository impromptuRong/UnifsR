function phylo = phyloseq(varargin)
%% Initialize Input Variables %%
import bioma.data.*
fieldlist = {'otu_table', 'phy_tree', 'sam_table', 'tax_table', 'Delimiter', 'seqdep'};
phylo = struct();
phylo.Delimiter = ',';
phylo.seqdep = true;
for i=1:2:length(varargin)
    field = varargin{i};
    if(ischar(field) && ismember(field, fieldlist))
        phylo.(field) = varargin{i+1};
    else
        error('phyloseq:InvalidInput', 'Invalid field Name');
    end
end

%% Read in Abundance Data Info %%
if(~isfield(phylo, 'otu_table'))
    error('phyloseq:argChk', 'Taxa Table is empty!');
else
    otu_table = phylo.otu_table;
    if(ischar(otu_table))
        otu_table = DataMatrix('File', otu_table, 'Delimiter', phylo.Delimiter);
    end
    rownames = regexprep(otu_table.RowNames, {'("|'')', '[^a-zA-Z0-9]'}, {'','_'});
    colnames = regexprep(otu_table.ColNames, {'("|'')', '[^a-zA-Z0-9]'}, {'','_'});
    otu_table = set(otu_table, 'RowNames', rownames, 'ColNames', colnames);
    if(length(phylo.seqdep)==1)
        sam_abu = phylo.seqdep * ones(size(otu_table,1),1);
    end
    if(islogical(phylo.seqdep))
        sam_abu = sum(otu_table, 2);
    end
end

%% Read in Phylogenetic Tree %%
if(~isfield(phylo, 'phy_tree'))
    error('phyloseq:argChk', 'Phylogenetic Tree is empty!');
else
    phy_tree = phylo.phy_tree;
    if(ischar(phy_tree))
        phy_tree = phytreeread(phy_tree);
    end
    phy_tree = reroot(phy_tree);
    nodename = get(phy_tree, 'NodeNames');
    nodename = regexprep(nodename, {'("|'')', '[^a-zA-Z0-9]'}, {'','_'});
    branch = get(phy_tree, 'Distance');
    branch(branch<=0) = 1e-8;
    phy_tree = phytree(get(phy_tree, 'Pointers'), branch, nodename);
end

%% Read in MetaData Info %%
if(~isfield(phylo, 'sam_table'))
    warning('phyloseq:argChk', 'Sample Table is empty!');
    sam_table.Sam_ID = rownames;
else
    sam_table = phylo.sam_table;
    if(ischar(sam_table))
        readinfo = csvimport(sam_table, 'delimiter', phylo.Delimiter);
        sam_table = struct;
        for i=1:size(readinfo,2)
            field = regexprep(char(readinfo(1,i)), {'("|'')', '[^a-zA-Z0-9]'}, {'','_'});
            tmp = readinfo(2:end,i);
            if ~iscellstr(tmp)
                tmp = cell2mat(tmp);
            else
                tmp = regexprep(tmp, '("|'')', '');
            end
            sam_table.(field) = tmp;
        end
    end
end
if(~isfield(sam_table, 'Sam_ID'))
    tmp = fieldnames(sam_table);
    sam_table.Sam_ID = sam_table.(tmp{1});
end

%% Read in Taxonomy Files %%
if(~isfield(phylo, 'tax_table'))
    warning('phyloseq:argChk', 'Taxa Table is empty!');
    tax_table.Tax_ID = colnames';
else
    tax_table = phylo.tax_table;
    if(ischar(tax_table))
        readinfo = csvimport(tax_table, 'delimiter', phylo.Delimiter);
        tax_table = struct;
        for i=1:size(readinfo,2)
            field = regexprep(char(readinfo(1,i)), {'("|'')', '[^a-zA-Z0-9]'}, {'','_'});
            tmp = readinfo(2:end,i);
            if ~iscellstr(tmp)
                tmp = cell2mat(tmp);
            else
                tmp = regexprep(tmp, {'("|'')', '[^a-zA-Z0-9]'}, {'','_'});
            end
            tax_table.(field) = tmp;
        end
    end
end
if(~isfield(tax_table, 'Tax_ID'))
    tmp = fieldnames(tax_table);
    tax_table.Tax_ID = tax_table.(tmp{1});
end
Tax_ID = tax_table.Tax_ID;
tax_table = structfun(@(x) [NaN;x], tax_table, 'UniformOutput', false);

%% Combine Info and Calculate Unifrac Structure %%
[common_tax, t_t, d_t] = intersect(get(phy_tree,'LeafNames'), otu_table.ColNames);
[~, i_t] = ismember(common_tax, Tax_ID);
[~, m_s, d_s] = intersect(sam_table.Sam_ID, otu_table.RowNames);
phylo.phy_tree = prune(phy_tree, setdiff(1:get(phy_tree,'NumLeaves'), t_t));
phylo.sam_table = structfun(@(x) x(m_s), sam_table, 'UniformOutput', false);
phylo.otu_table = otu_table(d_s, d_t);
sam_abu = sam_abu(d_s);
phylo.tax_table = structfun(@(x) x(i_t+1), tax_table, 'UniformOutput', false);
phylo.tax_table.Tax_ID = common_tax';
[~, phylo.edgebool, phylo.edgematrix] = fastUnifrac(phylo.otu_table, phylo.phy_tree, sam_abu);
[phylo.Nsample, phylo.Ntaxa] = size(phylo.otu_table);
phylo.Nexd = size(phylo.edgematrix.mat,2);
phylo = rmfield(phylo, 'Delimiter');
phylo = rmfield(phylo, 'seqdep');
% clearvars readinfo metaInfo dataInfo leafname branch tree rownames colnames;
% clearvars tmp i common_tax common_sam t_t d_t m_s d_s;
end
