library(mvtnorm)
library(fields)
p = 50
A = rmvnorm(n = p, sigma = diag(1,p))
x11();image.plot(A)
ACheatmap(A, center_value = 0, horizontal = F)
ACheatmap_hcl(A, horizontal = T)

n = 100
p = 6
data = rmvnorm(n = n, sigma = diag(1,p))
boxplot_from_matrix(data = data[,1:3], use_ggplot = T)

colnames(data) = c("uno","due","tre")
boxplot_from_matrix(data = data[,1:3], use_ggplot = F)

n1 = 20
n2 = 40
n3 = 10
v1 = rnorm(n = n1)
v2 = rnorm(n = n2, mean = 0.5)
v3 = rnorm(n = n3, sd = 0.75)
boxplot_from_vectors(v1,v2,v3, names = 1:3)

