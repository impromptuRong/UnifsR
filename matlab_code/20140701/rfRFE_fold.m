function [acc_tab, feature, model] = rfRFE_fold(cv_fold, ntree, mtry, extra_options, pattern)
kfold = length(cv_fold) - 1;
if isempty(ntree) || ntree < 1
    ntree = 500;
end
if isempty(mtry) || mtry < 1
    mtryf = @(x) floor(sqrt(size(x,2)));
else
    mtryf = @(x) mtry;
end

k = 1;
acck = zeros(kfold, 1);
parfor m = 1:kfold
    modk = classRF_train(cv_fold{m}{1}, cv_fold{m}{2}, ntree, mtryf(cv_fold{m}{1}), extra_options);
%     labelk = classRF_predict(cv_fold{m}{3}, modk);
%     confmat = confusionmat(cv_fold{m}{4}, labelk);
%     acck(m) = (confmat(1)+confmat(4))/sum(confmat(:)) * 100;
    acck(m) = 100 - modk.errtr(end,1) * 100;
end
disp(sprintf('Cross Validation Accuracy: %.2f%', mean(acck)));
cv(k) = mean(acck);

[CVacc, index] = max(cv);
mod = classRF_train(cv_fold{end}{1}, cv_fold{end}{2}, ntree, mtryf(cv_fold{end}{1}), extra_options);
pred = classRF_predict(cv_fold{end}{3}, mod);
confmat = confusionmat(cv_fold{end}{4}, pred);
TCacc = (confmat(1)+confmat(4))/sum(confmat(:)) * 100;
rankingCriteria = mod.importance;
nfea = length(rankingCriteria);
acc_tab = zeros(nfea,6);
feature = cell(nfea,1);
% model = cell(nfea,1);

k = 1;
acc_tab(k,:) = [length(rankingCriteria), CVacc(1), TCacc(1), ntree, mtryf(cv_fold{end}{1}), kfold];
feature{k} = 1:length(rankingCriteria);
% model{k} = mod;

while acc_tab(k,2) > max(acc_tab(:,2))-10 && length(rankingCriteria) > 2
    k = k + 1;
    sub_fea = rankingCriteria > quantile(rankingCriteria, 0.05);
    feature{k} = feature{k-1}(sub_fea);

    for m = 1:kfold+1
        cv_fold{m}{1} = cv_fold{m}{1}(:,sub_fea);
        cv_fold{m}{3} = cv_fold{m}{3}(:,sub_fea);
    end
    
    if ~strcmp(pattern, 'fast')
%         parfor k = 1:size(exp_grid,1)
%             par = exp_grid(k,:);
            k = 1;
            acck = zeros(kfold,1);
            parfor m = 1:kfold
                modk = classRF_train(cv_fold{m}{1}, cv_fold{m}{2}, ntree, mtryf(cv_fold{m}{1}), extra_options);
                labelk = classRF_predict(cv_fold{m}{3}, modk);
                confmat = confusionmat(cv_fold{m}{4}, labelk);
                acck(m) = (confmat(1)+confmat(4))/sum(confmat(:)) * 100;
            end
            disp(sprintf('Cross Validation Accuracy: %.2f%\n', mean(acck)));
            cv(k) = mean(acck);
%         end
        [CVacc, index] = max(cv);
    else
        acck = zeros(kfold,1);
        parfor m = 1:kfold
            modk = classRF_train(cv_fold{m}{1}, cv_fold{m}{2}, ntree, mtryf(cv_fold{m}{1}), extra_options);
%             labelk = classRF_predict(cv_fold{m}{3}, modk);
%             confmat = confusionmat(cv_fold{m}{4}, labelk);
%             acck(m) = (confmat(1)+confmat(4))/sum(confmat(:)) * 100;
            acck(m) = 100 - modk.errtr(end,1) * 100;
        end
        disp(sprintf('Cross Validation Accuracy: %.2f%', mean(acck)));
        CVacc = mean(acck);
    end
    mod = classRF_train(cv_fold{end}{1}, cv_fold{end}{2}, ntree, mtryf(cv_fold{end}{1}), extra_options);
    pred = classRF_predict(cv_fold{end}{3}, mod);
    confmat = confusionmat(cv_fold{end}{4}, pred);
    TCacc = (confmat(1)+confmat(4))/sum(confmat(:)) * 100;
    rankingCriteria = mod.importance;
    acc_tab(k,:) = [length(rankingCriteria), CVacc(1), TCacc(1), ntree, mtryf(cv_fold{end}{1}), kfold];
%    model{k} = mod;
end
acc_tab = acc_tab(1:k,:);
feature = feature(1:k);
model = model(1:k);

end

