#fast FET test
#generate 1000 random contingency tables
a = round(sample(0:1000, 1000))
b = round(sample(0:5000, 1000))
c = round(sample(0:5000, 1000))
d = round(sample(0:20000, 1000))


test = fastFET(a,b,c,d, alternative = "greater")
true = numeric(length(a))
for(i in 1:length(a)) true[i] = fisher.test(matrix(c(a[i],c[i],b[i],d[i]),ncol = 2),alternative = 'greater')$p.value

print(sum(abs(true-test)))
