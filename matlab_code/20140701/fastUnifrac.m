function [edge_list, edge_bool, edge_matrix, dist_list] = fastUnifrac(data, tree, seqdep, method)
%% Check Input error %%
import bioma.data.*
if(nargin<4)||isempty(method),      method='Edge';      end
if(nargin<3)||isempty(seqdep),      seqdep=sum(data,2); end
method = ismember({'Edge','dPCoA','non_w','w.non_nor','w.nor'}, method);
if(~sum(method))
    warning('fastUni:methodCk', 'Methods are not Identified!\nOnly return edgeMatrix and edgeList.\nProper Methods are: Edge; dPCoA; non_w; w.non_nor; w.nor');
end

%% Get info anc check structures of data and tree %%
edgedist = get(tree, 'Distance');
edgedist(edgedist<=0) = 1e-4;
if(isempty(edgedist)),      error('fastUni:treeDisCk', 'Tree has no branch lengths!');  end
Ntaxa = get(tree, 'NumLeaves');
Nnode = get(tree, 'NumNodes');
Nedge = get(tree, 'NumBranches');
if((Ntaxa+Nedge)~=Nnode),   error('fastUni:treeStrCk', 'Tree structure error!');        end
leaf = get(tree, 'LeafNames');
node = get(tree,'NodeNames');
branch = get(tree, 'Pointers');
if(~all(ismember(leaf,data.ColNames))), error('fastUni:dataTreMt', 'Data miss Tree leaves!');   end
data = data./repmat(seqdep, 1, Ntaxa);
data = data(:, leaf);

%% create edge_list and edge_matrix %%
%    function sons = extendedge(edge)
%        if(edge <= Ntaxa)
%            sons = edge;
%        else
%            edge = edge - Ntaxa;
%            sons = branch(edge, :);
%            while(any(sons>Ntaxa))
%                sons = [reshape(branch((sons(sons>Ntaxa)-Ntaxa),:),1,[]), sons(sons<=Ntaxa)];
%            end
%        end
%    end
%    edge_list = arrayfun(@extendedge, 1:Nnode, 'UniformOutput', false);

parfor edgex = 1:Nnode
    if(edgex <= Ntaxa)
        sons = edgex;
    else
        edge = edgex - Ntaxa;
        sons = branch(edge, :);
        while(any(sons>Ntaxa))
            sons = [reshape(branch((sons(sons>Ntaxa)-Ntaxa),:),1,[]), sons(sons<=Ntaxa)];
        end
    end
    edge_list{edgex} = sons;
end

edge_cell = cellfun(@(x) ismember(1:Ntaxa, x)', edge_list, 'UniformOutput', false);
% edge_bool = DataMatrix(cell2mat(edge_cell), leaf, node);
edge_bool.mat = sparse(cell2mat(edge_cell));
edge_bool.rownames = leaf;
edge_bool.colnames = node;
% edge_matrix = DataMatrix(double(data)*double(edge_bool), data.RowNames, node);
edge_matrix.mat = sparse(double(data))*edge_bool.mat;
edge_matrix.rownames = data.RowNames;
edge_matrix.colnames = node;

%% calculate distance matrix %%
dist_list = struct;
if(method(2))
    D = pdist(edge_matrix.mat, @(x,J) sqrt(double((repmat(x,size(J,1),1)-J)).^2*edgedist));
    dist_list.dPCoA = DataMatrix(squareform(D), edge_matrix.rownames, edge_matrix.rownames);
end
if(method(3))
    D = pdist(edge_matrix.mat, @(x,J) (abs(double(repmat(x>0,size(J,1),1)-(J>0)))*edgedist)./(abs(double(repmat(x>0,size(J,1),1)|(J>0)))*edgedist));
    dist_list.non_w = DataMatrix(squareform(D), edge_matrix.rownames, edge_matrix.rownames);
end
if(method(4))
    D = pdist(edge_matrix.mat, @(x,J) abs(double(repmat(x,size(J,1),1)-J))*edgedist);
    dist_list.w.non_nor = DataMatrix(squareform(D), edge_matrix.rownames, edge_matrix.rownames);
end
if(method(5))
    D = pdist(edge_matrix.mat, @(x,J) (abs(double(repmat(x,size(J,1),1)-J))*edgedist)./(abs(double(repmat(x,size(J,1),1)+J))*edgedist));
    dist_list.w.nor = DataMatrix(squareform(D), edge_matrix.rownames, edge_matrix.rownames);
end
end
