library(MASS)



ii <- c(2, 3, 4, 3, 4, 4)
jj <- c(1, 1, 1, 2, 2, 3)
weights <- 1:6
delta <- 1:6

ei <- function(i, n) {
  return(ifelse(i == (1:n), 1, 0))
}

aij <- function(i, j, n) {
  df <- ei(i, n) - ei(j, n)
  return(outer(df, df))
}

makeU <- function(n) {
  ij <- 1
  m <- n * (n - 1) / 2
  u <- matrix(0, m, m) 
  for (j in 1:(n - 1)) {
    for (i in (j + 1):n) {
      a <- aij(i, j, n)
      kl <- 1
      for (l in 1:(n - 1)) {
        for (k in (l + 1):n) {
          b <- aij(k, l, n)
          u[ij, kl] <- sum(a * b)
          kl <- kl + 1
        }
      }
      ij <- ij + 1
    }
  }
  return(u)
}

triangular2SDC <- function(x) {
  k <- 1
  n <- makeN(length(x))
  s <- matrix(0, n, n)
  for (j in 1:(n - 1)) {
    for(i in (j + 1):n) {
      s[i, j] <- s[j, i] <- -x[k]
      k <- k + 1
    }
  }
  diag(s) <- -rowSums(s)
  return(s)
}

triangular2Hollow <- function(x) {
  k <- 1
  n <- makeN(length(x))
  s <- matrix(0, n, n)
  for (j in 1:(n - 1)) {
    for(i in (j + 1):n) {
      s[i, j] <- s[j, i] <- x[k]
      k <- k + 1
    }
  }
  return(s)
}  

doubleCenter <- function(x) {
  rs <- apply(x, 1, mean)
  ss <- mean(x)
  return(x - outer(rs, rs, "+") + ss)
}

# mPrint() formats a matrix (or vector, or scalar) of numbers
# for printing 

mPrint <- function(x,
                   digits = 6,
                   width = 8,
                   format = "f",
                   flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}

makeN <- function(m) {
  return(as.integer((1 + sqrt(1 + 8 * m)) / 2))
}

torgerson <- function(delta, weights, p = 2) {
  n <- makeN(length(delta))
  u <- makeU(n)
  th <- lsfit(u, delta ^ 2, wt = weights, intercept = FALSE)$coef
  mPrint(triangular2SDC(th))
  mPrint(triangular2SDC(solve(u, delta)))
}

makeC <- function(weights) {
  k <- 1
  m <- length(weights)
  n <- makeN(m)
  s <- matrix(0, n ^ 2, n ^ 2)
  t <- matrix(0, n ^ 2, m)
  k <- 1
  for (j in 1:(n - 1)) {
    for (i in (j + 1):n) {
      t[, k] <- as.vector(aij(i, j, n))
      k <- k + 1
    }
  }
  u <- crossprod(t)
  u <- u * outer(sqrt(weights), sqrt(weights))
  return(list(s = s, u = u))
}

perronRoot <- function(b,
                       lbd = 0.0,
                       itmax = 100,
                       eps = 1e-10,
                       verbose = TRUE) {
  n <- nrow(b)
  a  <- b + lbd * diag(n)
  itel = 1
  
  repeat {
    r <- rowSums(a)
    rmax <- max(r)
    rmin <- min(r)
    if (verbose) {
      cat("itel", formatC(itel, digits = 3, format = "d"),
          "rmin", formatC(rmin, digits = 10, format = "f"),
          "rmax", formatC(rmax, digits = 10, format = "f"),
          "\n")
    }
    if (((rmax - rmin) < eps) || (itel == itmax)) {
      break
    }
    a <- a * t(outer(r, r, "/"))
    itel <- itel + 1
  }
  return(((rmin + rmax) / 2.0) - lbd)
}
