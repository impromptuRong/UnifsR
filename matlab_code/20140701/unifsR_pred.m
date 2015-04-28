function [result, score] = unifsR_pred(model, TEST, gclass)
score = [ones(size(TEST,1),1), double(TEST)]*model.w+10^-20;
pred = nominal([-1;1;sign(score);sign(-score)], model.label);
pred = pred((1:size(score,1))+2);

switch gclass
    case 'nominal'
        result = pred;
    case 'char'
        result = char(pred);
    case 'cell'
        result = cellstr(pred);
    otherwise
        result = sign(score);
end