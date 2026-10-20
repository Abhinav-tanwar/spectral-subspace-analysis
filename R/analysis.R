null_basis <- function(A, tol = 1e-10) {
    q <- qr(t(A), tol = tol, LAPACK = FALSE)
    Q <- qr.Q(q, complete = TRUE)
    if (q$rank == ncol(A)) 
        return(matrix(numeric(), ncol(A), 0))
    Q[, seq.int(q$rank + 1, ncol(A)), drop = FALSE]
}

run_analysis <- function(root, out) {
    set.seed(3)
    square <- matrix(sample(-10:10, 100, TRUE), 10, byrow = TRUE)
    se <- eigen(square)
    inv <- solve(square)
    square_rank <- qr(square)$rank
    image_basis <- qr.Q(qr(square))[, seq_len(square_rank), drop = FALSE]
    set.seed(6)
    rectangular <- matrix(sample(-10:10, 32, TRUE), 4, 8)
    rectangular_null <- null_basis(rectangular)
    set.seed(7)
    v <- runif(5)
    set.seed(8)
    u <- runif(3)
    outer_product <- v %o% u
    symmetric <- crossprod(square)
    sy <- eigen(symmetric, symmetric = TRUE)
    write.csv(image_basis, file.path(out, "image-basis.csv"), row.names = FALSE)
    write.csv(rectangular_null, file.path(out, "rectangular-null-basis.csv"), row.names = FALSE)
    set.seed(530302)
    A <- matrix(rnorm(15), 3, 5)
    b <- c(2, -1, 4)
    N <- null_basis(A)
    xp <- MASS::ginv(A) %*% b
    z <- c(0.3, -0.7)
    B <- matrix(c(0, -1, 1, 0), 2, byrow = TRUE)
    e <- eigen(B)
    reconstructed <- e$vectors %*% diag(e$values) %*% solve(e$vectors)
    set.seed(9)
    large <- matrix(runif(128^2), 128)
    values <- eigen(large, only.values = TRUE)$values
    write.csv(data.frame(real = Re(values), imaginary = Im(values)), file.path(out, "eigenvalues.csv"), 
        row.names = FALSE)
    list(square_rank = square_rank, square_nullity = ncol(null_basis(square)), determinant_inverse_relative_error = abs(det(inv) * 
        det(square) - 1), image_orthonormal_error = max(abs(crossprod(image_basis) - diag(square_rank))), 
        rank_transpose_agrees = qr(rectangular)$rank == qr(t(rectangular))$rank, rectangular_rank_nullity = qr(rectangular)$rank + 
            ncol(rectangular_null), outer_product_rank = qr(outer_product)$rank, symmetric_orthonormal_error = max(abs(crossprod(sy$vectors) - 
            diag(10))), square_eigen_residual = max(Mod(square %*% se$vectors - se$vectors %*% 
            diag(se$values))), rank = qr(A)$rank, nullity = ncol(N), rank_nullity = qr(A)$rank + 
            ncol(N), null_residual = max(abs(A %*% N)), solution_family_residual = max(abs(A %*% 
            (xp + N %*% z) - b)), eigen_residual = max(Mod(B %*% e$vectors - e$vectors %*% 
            diag(e$values))), reconstruction_residual = max(Mod(B - reconstructed)), largest_eigenvalue_magnitude = max(Mod(values)), 
        smallest_eigenvalue_magnitude = min(Mod(values)))
}

