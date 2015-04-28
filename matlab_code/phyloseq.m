function phylo = phyloseq(varargin)
%% Initialize Input Variables %%
import bioma.data.*
fieldlist = {'otu_table', 'phy_tree', 'sam_table', 'tax_table', 'Delimiter', 'seqdep', 'trim', 'unifrac'};
phylo = struct();
phylo.Delimiter = ',';
phylo.seqdep = true;
phylo.unifrac = true;
phylo.trim = true;
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
    if(isempty(otu_table.RowNames))
        otu_table = set(otu_table(:,2:end), 'RowNames', cellstr(num2str(double(otu_table(:, 1)))));
    end
    rownames = regexprep(otu_table.RowNames, {'("|'')', '[^a-zA-Z0-9]'}, {'','_'});
    colnames = regexprep(otu_table.ColNames, {'("|'')', '[^a-zA-Z0-9]'}, {'','_'});
    otu_table = double(otu_table);
    if(islogical(phylo.seqdep))
        if(phylo.seqdep)
            sam_abu = sum(otu_table, 2);
        else
            sam_abu = ones(size(otu_table,1),1);
        end
    elseif(length(phylo.seqdep)==1)
        sam_abu = phylo.seqdep * ones(size(otu_table,1),1);
    else
        sam_abu = phylo.seqdep;
    end
end

%% Trim otu_table based on cut-off %%
if(islogical(phylo.trim) && phylo.trim)
    phylo.trim = [0; 0];
end
if(~islogical(phylo.trim))
    rowIndex = mean(abs(otu_table),2) > phylo.trim(1);
    colIndex = mean(abs(otu_table),1) > phylo.trim(2);
    otu_table = otu_table(rowIndex, colIndex);
    rownames = rownames(rowIndex);
    colnames = colnames(colIndex);
    sam_abu = sam_abu(rowIndex);
end

%% Read in MetaData Info %%
if(~isfield(phylo, 'sam_table'))
    warning('Sample Table is empty!');
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
if(isnumeric(sam_table.Sam_ID))
    sam_table.Sam_ID = cellstr(num2str(sam_table.Sam_ID));
end

%% Sync sam_table and otu_table %%
[rownames, m_s, d_s] = intersect(sam_table.Sam_ID, rownames);
sam_table = structfun(@(x) x(m_s), sam_table, 'UniformOutput', false);
otu_table = otu_table(d_s, :);
sam_abu = sam_abu(d_s);

%% Read in Phylogenetic Tree %%
if(~isfield(phylo, 'phy_tree'))
    warning('Phylogenetic Tree is empty!');
    phylo.unifrac = false;
else
    phy_tree = phylo.phy_tree;
    if(ischar(phy_tree))
        phy_tree = phytreeread(phy_tree);
    end
    phy_tree = reroot(phy_tree);
    nodename = get(phy_tree, 'NodeNames');
    nodename = regexprep(nodename, {'("|'')', '[^a-zA-Z0-9]'}, {'','_'});
    branch = get(phy_tree, 'Distance');
    branch(branch<0) = 0;
    phy_tree = phytree(get(phy_tree, 'Pointers'), branch, nodename);
end

%% Sync phy_tree and otu_table %%
if(isfield(phylo, 'phy_tree'))
    [colnames, t_t, d_t] = intersect(get(phy_tree,'LeafNames'), colnames);
    phy_tree = prune(phy_tree, setdiff(1:get(phy_tree,'NumLeaves'), t_t));
    otu_table = otu_table(:, d_t);
end

%% Read in Taxonomy Files %%
if(~isfield(phylo, 'tax_table'))
    warning('Taxa Table is empty!');
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
if(isnumeric(tax_table.Tax_ID))
    tax_table.Tax_ID = cellstr(num2str(tax_table.Tax_ID));
end
Tax_ID = tax_table.Tax_ID;
tax_table = structfun(@(x) [NaN;x], tax_table, 'UniformOutput', false);

%% Sync tax_table with otu_table %%
[~, i_t] = ismember(colnames, Tax_ID);
tax_table = structfun(@(x) x(i_t+1), tax_table, 'UniformOutput', false);

%% Calculate Unifrac %%
phylo.otu_table = DataMatrix(otu_table, 'RowNames', rownames, 'ColNames', colnames);
phylo.sam_table = sam_table;
phylo.tax_table = tax_table;
phylo.phy_tree = phy_tree;

if(phylo.unifrac)
    [~, phylo.edgebool, phylo.edgematrix] = fastUnifrac(phylo.otu_table, phylo.phy_tree, sam_abu);
    phylo.Nexd = size(phylo.edgematrix.mat,2);
end
[phylo.Nsample, phylo.Ntaxa] = size(phylo.otu_table);
phylo.seqdep = sam_abu;
phylo = rmfield(phylo, {'Delimiter', 'unifrac'});
% clearvars readinfo metaInfo dataInfo leafname branch tree rownames colnames;
% clearvars tmp i common_tax common_sam t_t d_t m_s d_s;
end
