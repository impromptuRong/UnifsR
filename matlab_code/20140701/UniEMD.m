matlabpool local 8;
addpath(genpath(pwd));
import bioma.data.*

%% Read in Pairwise DNA distance Matrix %%
load('oralotu_phylo.mat');
dl = struct;
dl.samples = get(phylo.data, 'RowNames');
dl.species = get(phylo.data, 'ColNames');
dl.meta = phylo.meta;

distDNA = DataMatrix('File', '0.dist.csv', 'Delimiter', ',');
otuName = regexprep(get(distDNA, 'RowNames'), '"', '');
otuName = regexprep(otuName, '[^a-zA-Z0-9]', '_');
distDNA = set(distDNA, 'RowNames', otuName, 'ColNames', otuName);
distDNA = distDNA(dl.species, dl.species);

%% Calculate Upper level distance Matrix %%
taxaLevel = {'Domain','Phylum','Class','Order','Family','Genus','OTU'};
parIndex = combnk(1:phylo.Nsample, 2);

for k=1:length(taxaLevel)
    level = taxaLevel{k};
    D = zeros(1,size(parIndex,1));
    taxaIndex = cellstr(getlevels(nominal(phylo.taxa.(level))));
    for i=1:length(taxaIndex)
        otuIndex = strcmp(taxaIndex(i),phylo.taxa.(level));
        abuMat = double(phylo.data(dl.samples, otuIndex));
        dnaMat = double(distDNA(otuIndex, otuIndex));

        tmp = zeros(1,size(parIndex,1));
        parfor j=1:size(parIndex,1)
            tmp(j) = emd_hat_gd_metric_mex(abuMat(parIndex(j,1),:)', abuMat(parIndex(j,2),:)', dnaMat, 0);
        end
        D = D + tmp;
    end
    dl.EMD.(level) = squareform(D);
    dmwrite(DataMatrix(dl.EMD.(level), dl.samples, dl.samples), strcat(level, '.EMD.csv'), 'Delimiter', ',');
%    DataMatrix(squareform(D), dl.samples, dl.samples);
end

%% Calculate Unifrac distance Matrix %%
load('oralotu_phylo.mat');
edgedist = get(phylo.tree, 'Distance');
edgedist(edgedist<=0) = 1e-4;
D = pdist(phylo.edgematrix.mat, @(x,J) double(abs(repmat(x>0,size(J,1),1)-(J>0)))*edgedist);
uni.OTU.non_w = squareform(D);

D = pdist(phylo.edgematrix.mat, @(x,J) double(abs(repmat(x,size(J,1),1)-J))*edgedist);
uni.OTU.w_non_nor = squareform(D);

D = pdist(phylo.edgematrix.mat, @(x,J) (double(abs(repmat(x,size(J,1),1)-J))*edgedist)./(double(abs(repmat(x,size(J,1),1)+J))*edgedist));
uni.OTU.w_nor = squareform(D);

dmwrite(DataMatrix(dl.Unifrac.OTU.non_w, dl.samples, dl.samples), 'OTU.Unifrac.non_w.csv', 'Delimiter', ',');
dmwrite(DataMatrix(dl.Unifrac.OTU.w_non_nor, dl.samples, dl.samples), 'OTU.Unifrac.w_non_nor.csv', 'Delimiter', ',');
dmwrite(DataMatrix(dl.Unifrac.OTU.w_nor, dl.samples, dl.samples), 'OTU.Unifrac.w_nor.csv', 'Delimiter', ',');

matlabpool close;

