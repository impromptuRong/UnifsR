function dist = unifsR_dis(x, y, beta1, beta2, alpha, npar)
p = full(abs(x-y));
if(alpha == 0)
    dist = sqrt(p.^2 * beta2 / npar);
    return;
end
if(alpha == 1)
    dist = max(p'.*beta1);
    return;
end
noZero = beta1~=0;
p = p(noZero);
beta1 = beta1(noZero);
beta2 = beta2(noZero);

[Rseq, index] = sort([p'.*beta1;0], 'descend');
R = Rseq(1);  
M = 0;
i = 1;
while((1/alpha-1)^2*R^2*npar > M)
    i = i + 1;
%     while(Rseq(i) > R-10^-4)
%         i = i + 1;
%     end
    R = Rseq(i);
    M = sum((Rseq(1:i)-R).^2.*beta2(index(1:i))./(beta1(index(1:i)).^2));
end
list = index(1:i-1);
p = p(list);
beta1 = beta1(list);
beta2 = beta2(list);
a = sum(beta2./(beta1.^2)) - npar*(1/alpha-1)^2;
b = p*(beta2./beta1);
c = p.^2*beta2;
dist = full((b-sqrt(b^2-a*c))/a/alpha);
end