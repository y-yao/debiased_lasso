function reject = TwoSampleTTest(sample1, sample2, alpha_value)
% two-tailed t-test
% variances are assumed to be unequal

mean1 = mean(sample1); mean2 = mean(sample2);

var1 = var(sample1); var2 = var(sample2);

n1 = length(sample1); n2 = length(sample2);

t = (mean1 - mean2) / sqrt(var1^2/n1 + var2^2/n2);

if t > tinv(1 - alpha_value/2, min(n1, n2))
    reject = 1;
else
    reject = 0;
end
end