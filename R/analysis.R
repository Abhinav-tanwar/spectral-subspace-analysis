null_basis <- function(A,tol=1e-10) {
 q <- qr(t(A),tol=tol,LAPACK=FALSE); Q <- qr.Q(q,complete=TRUE)
 if(q$rank == ncol(A)) return(matrix(numeric(),ncol(A),0))
 Q[,seq.int(q$rank+1,ncol(A)),drop=FALSE]
}
run_analysis <- function(root,out) {
 set.seed(530302); A <- matrix(rnorm(15),3,5); b <- c(2,-1,4)
 N <- null_basis(A); xp <- MASS::ginv(A)%*%b; z <- c(.3,-.7)
 B <- matrix(c(0,-1,1,0),2,byrow=TRUE); e <- eigen(B)
 reconstructed <- e$vectors%*%diag(e$values)%*%solve(e$vectors)
 set.seed(9); large <- matrix(runif(128^2),128)
 values <- eigen(large,only.values=TRUE)$values
 write.csv(data.frame(real=Re(values),imaginary=Im(values)),file.path(out,"eigenvalues.csv"),row.names=FALSE)
 list(rank=qr(A)$rank,nullity=ncol(N),rank_nullity=qr(A)$rank+ncol(N),null_residual=max(abs(A%*%N)),solution_family_residual=max(abs(A%*%(xp+N%*%z)-b)),eigen_residual=max(Mod(B%*%e$vectors-e$vectors%*%diag(e$values))),reconstruction_residual=max(Mod(B-reconstructed)),largest_eigenvalue_magnitude=max(Mod(values)),smallest_eigenvalue_magnitude=min(Mod(values)))
}
