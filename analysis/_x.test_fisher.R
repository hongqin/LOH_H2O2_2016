mat = matrix(nrow=2,ncol=2)
mat[1,1]=30
mat[2,2]=45
mat[1,2]=5
mat[2,1]=6
mat;
fisher.test(mat)

mat = mat/sum(mat)
mat;
fisher.test(mat)

#conclusion. R's defaulty fisher.test does not handle proportion input. 
