function [acc_tab, feature, model] = svmRFE(trainy, trainX, cmd, parlist, pattern)
trainX = double(trainX);
kernel = str2double(cmd(regexp(cmd, '-t ')+3));
if isempty(pattern)
    pattern = 'fast';
end
if ~isfield(parlist,'c') || isempty(parlist.c)
    parlist.c = 2.^(-5:15);
end
if ~isfield(parlist,'g') || isempty(parlist.g)
    if kernel~=0
        parlist.g = 2.^(-15:3);
    else
        parlist.g = 1;
    end
end
if ~isfield(parlist,'d') || isempty(parlist.d)
    if kernel==1
        parlist.d = 2:5;
    else
        parlist.d = 1;
    end
end
[g1, g2, g3] = meshgrid(parlist.c, parlist.g, parlist.d);
exp_grid = [g1(:), g2(:), g3(:)];

parfor k = 1:size(exp_grid,1)
    par = exp_grid(k,:);
    cv(k) = libsvmtrain(trainy, trainX, ['-s 0 -v 10 -q ', cmd, ' -c ', num2str(par(1)), ' -g ', num2str(par(2)), ' -d ', num2str(par(3))]);
end
% fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', exp_grid(index,1), exp_grid(index,2), cv(index), bestc, bestg, bestcv);
[CVacc, index] = max(cv);
cmd0 = ['-s 0 -v 10 -q ', cmd, ' -c ', num2str(exp_grid(index,1)), ' -g ', num2str(exp_grid(index,2)), ' -d ', num2str(exp_grid(index,3))];
mod = libsvmtrain(trainy, trainX, regexprep(cmd0, '-v (\d+)', ''));
[~, TCacc, ~] = libsvmpredict(trainy, trainX, mod, '-q');
rankingCriteria = svmWeights(mod);
nfea = length(rankingCriteria);
acc_tab = zeros(nfea,6);
feature = cell(nfea,1);
model = cell(nfea,1);

k = 1;
acc_tab(k,:) = [length(rankingCriteria), CVacc, TCacc(1), exp_grid(index,:)];
feature{k} = 1:length(rankingCriteria);
model{k} = mod;

while acc_tab(k,2) > max(acc_tab(:,2))-10 && length(rankingCriteria) > 2
    k = k + 1;
    sub_fea = rankingCriteria > quantile(rankingCriteria,0.05);
    feature{k} = feature{k-1}(sub_fea);
    trainX = trainX(:,sub_fea);
    
    if ~strcmp(pattern, 'fast')
        parfor k = 1:size(exp_grid,1)
            par = exp_grid(k,:);
            cv(k) = libsvmtrain(trainy, trainX, ['-s 0 -v 10 -q ', cmd, ' -c ', num2str(par(1)), ' -g ', num2str(par(2)), ' -d ', num2str(par(3))]);
        end
        [CVacc, index] = max(cv);
        cmd0 = ['-s 0 -v 10 -q ', cmd, ' -c ', num2str(exp_grid(index,1)), ' -g ', num2str(exp_grid(index,2)), ' -d ', num2str(exp_grid(index,3))];
    else
        CVacc = libsvmtrain(trainy, trainX, cmd0);
    end
    mod = libsvmtrain(trainy, trainX, regexprep(cmd0, '-v (\d+)', ''));
    [~, TCacc, ~] = libsvmpredict(trainy, trainX, mod, '-q');
    rankingCriteria = svmWeights(mod);
    acc_tab(k,:) = [length(rankingCriteria), CVacc, TCacc(1), exp_grid(index,:)];
    model{k} = mod;
end

acc_tab = acc_tab(1:k,:);
feature = feature(1:k);
model = model(1:k);

end


function rankingCriteria = svmWeights(model)
%% weights and Criteria of the hiperplane
kernel = model.Parameters(2);
nfea = size(model.SVs, 2);
SV = model.SVs;
coef = model.sv_coef;
if kernel==0                % linear kernel: u'*v
    rankingCriteria = (coef' * SV).^2;
elseif kernel==1
    d = model.Parameters(3);
    gamma = model.Parameters(4);
elseif kernel==2            % RBF kernel: exp(-gamma*|u-v|^2)
    gamma = model.Parameters(4);
    K0 = pdist(SV).^2;
    w0 = coef' * exp(squareform(-gamma*K0)) * coef;
    parfor f = 1:nfea
        Kf = K0 - pdist(SV(:,f)).^2;
        rankingCriteria(f) = (w0 - coef' * exp(squareform(-gamma*Kf)) * coef)/2;
    end
end
end


% 
% 	else{
% 		start <- c(1, cumsum(model$nSV)+1)
% 		start <- start[-length(start)]
% 
% 		W <- matrix(0,ncol(model$SV),choose(model$nclasses,2))
% 		count <- 1
% 		for(i in 1:(model$nclasses-1)){
% 			for(j in (i+1):model$nclasses){
% 				## ranges for class i and j:
% 				ri <- start[i]:(start[i] + model$nSV[i] - 1)
% 				rj <- start[j]:(start[j] + model$nSV[j] - 1)
% 				## coefs and SV for (i,j):
% 				coefs <- c(model$coefs[ri, j-1], model$coefs[rj, i])
% 				SV <- data.matrix(model$SV[c(ri,rj),])
% 				if(model$kernel==0){
% 					w <- t(coefs) %*% SV
% 					W[,count] <- w * w
% 				}
% 				if(model$kernel==1){
% 				}
% 				if(model$kernel==2){
% 					for(nf in 1:ncol(model$SV)){
% 						KMat <- (coefs %*% t(coefs)) * kernelMatrix(rbfdot(sigma=model$gamma),SV[,-nf],SV[,-nf])
% 						W[nf,count] <- -sum(KMat)
% 					}
% 				}
% 				count <- count+1
% 			}
% 		}
% 		rankingCriteria <- rowMeans(W)
% 	}
% 	return(rankingCriteria)
% }
