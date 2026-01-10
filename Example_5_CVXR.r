## R script using CVXR to generate designs and Figure 1 from the paper 
## Optimal Exact Designs of Multiresponse Experiments under Linear and Sparsity Constraints (https://arxiv.org/abs/2507.04713).
## The choice of solver can be adjusted by the argument `solver` in the function CRLAS
## Available options are "ECOS_BB", "GUROBI" (if installed) and "MOSEK" (if installed)
## We strongly recommend using either "GUROBI" or "MOSEK"

library(OptimalDesign)

#--------------------------------------------------------
#---------------- helper functions ----------------------
#--------------------------------------------------------

MISOCP_m4 <- function(Fx, b1 = numeric(0), A1 = NULL, b2 = numeric(0), 
          A2 = NULL, b3 = numeric(0), A3 = NULL, bin = FALSE, type = c("exact", "approximate"), 
          gap = NULL, t.max = Inf, solver, mi_max_iters = 3000, 
          mi_abs_eps = 1e-5,  mi_rel_eps = 1e-5, verbose = FALSE) {
  
  type <- match.arg(type)
  Fx <- as.matrix(Fx)
  n <- nrow(Fx)
  m <- ncol(Fx)
  if (m != 4) stop("This function is specialized to m = 4 (ncol(Fx) must be 4).")
  
  if (!requireNamespace("CVXR", quietly = TRUE)) stop("Please install.packages('CVXR').")
  if (!requireNamespace("ECOSolveR", quietly = TRUE)) stop("Please install.packages('ECOSolveR').")
  library(CVXR)
  
  # --- decision variables ---
  if (type == "exact") {
    w <- if (isTRUE(bin)) Variable(n, boolean = TRUE) else Variable(n, integer = TRUE)
  } else {
    w <- Variable(n)
  }
  
  # Lift variables as in od_D_MISOCP / MISOCP_gurobi
  U <- Variable(n, 4)      # >= 0
  V <- Variable(n, 4)      # free
  R <- Variable(4, 4)      # will be constrained upper-triangular
  
  # Objective scalar: proportional to (prod diag(R))^(1/4)
  tau <- Variable(1)
  
  cons <- list(w >= 0, U >= 0, tau >= 0)
  
  # For approximate+bin, mimic the "ub=1" behavior
  if (type == "approximate" && isTRUE(bin)) cons <- c(cons, list(w <= 1))
  
  # Linear side constraints A1 w <= b1, A2 w >= b2, A3 w = b3
  if (length(b1) > 0) {
    if (is.null(A1)) stop("b1 provided but A1 is NULL.")
    cons <- c(cons, list(as.matrix(A1) %*% w <= b1))
  }
  if (length(b2) > 0) {
    if (is.null(A2)) stop("b2 provided but A2 is NULL.")
    cons <- c(cons, list(as.matrix(A2) %*% w >= b2))
  }
  if (length(b3) > 0) {
    if (is.null(A3)) stop("b3 provided but A3 is NULL.")
    cons <- c(cons, list(as.matrix(A3) %*% w == b3))
  }
  
  # --- structural constraints (m=4), mirroring MISOCP_gurobi ---
  # R = 0.5 * t(V) %*% Fx
  cons <- c(cons, list(R == 0.5 * t(V) %*% Fx))
  
  # R upper triangular: zero below diagonal
  for (col in 1:3) for (row in (col+1):4) cons <- c(cons, list(R[row, col] == 0))
  
  # Diagonal linkage: R[i,i] = sum_j U[j,i]
  for (i in 1:4) cons <- c(cons, list(R[i, i] == sum_entries(U[, i])))
  
  # SOC product constraints for each (j,i):
  # || (U[j,i] - w[j], V[j,i]) ||_2 <= U[j,i] + w[j]
  # which implies V[j,i]^2 <= 4 * U[j,i] * w[j] when U,w >= 0.
  for (i in 1:4) {
    for (j in 1:n) {
      cons <- c(cons, list(
        norm2(vstack(U[j, i] - w[j], V[j, i])) <= (U[j, i] + w[j])
      ))
    }
  }
  
  # --- geometric mean of the 4 diagonal entries of R via 3 SOC constraints ---
  d1 <- R[1, 1]; d2 <- R[2, 2]; d3 <- R[3, 3]; d4 <- R[4, 4]
  g12 <- Variable(1); g34 <- Variable(1)
  cons <- c(cons, list(g12 >= 0, g34 >= 0))
  
  # g12 <= sqrt(d1*d2)
  cons <- c(cons, list(norm2(vstack(d1 - d2, 2 * g12)) <= d1 + d2))
  # g34 <= sqrt(d3*d4)
  cons <- c(cons, list(norm2(vstack(d3 - d4, 2 * g34)) <= d3 + d4))
  # tau <= sqrt(g12*g34) = (d1*d2*d3*d4)^(1/4)
  cons <- c(cons, list(norm2(vstack(g12 - g34, 2 * tau)) <= g12 + g34))
  
  prob <- Problem(Maximize(tau), cons)
  
  if (type == "exact") {
    # ECOS_BB solver options are exposed in CVXR (mi_max_iters, mi_abs_eps, mi_rel_eps)
    res <- solve(prob, solver = solver, mi_max_iters = mi_max_iters,
                 mi_abs_eps = mi_abs_eps, mi_rel_eps = mi_rel_eps, verbose = verbose)
  } else {
    res <- solve(prob, solver = solver, verbose = verbose)
  }
  
  w_best <- res$getValue(w)
  if (type == "exact" && is.numeric(w_best)) w_best <- round(w_best)
  
  list(w.best = as.numeric(w_best), status = res$status, Phi.best = res$value)
}


ContF <- function(x, a, b){
  # elementary information matrices for continuation ratio (CR) model in dose x
  # x: the dose
  # a,b: parameter values of the CR model a=(a1,a2), b=(b1,b2)
  
  e1 <- exp(a[1]+b[1]*x)
  e2 <- exp(a[2]+b[2]*x)
  
  c1 <- e2 / ((1+e2)^2*(1+e1))
  c2 <- e1 / (1+e1)^2
  
  M <- rbind(c(1,x,0,0), c(x,x^2,0,0), c(0,0,1,x), c(0,0,x,x^2))
  M[1:2, ] <- M[1:2, ] * c1
  M[3:4, ] <- M[3:4, ] * c2
  
  E <- eigen(M)
  
  Fx <- matrix(0, nrow=2, ncol=4)
  for (i in 1:2) Fx[i, ] <- sqrt(E$values[i]) * E$vectors[ ,i]
  
  Fx
  
}

resp <- function(x, a, b){
  # the response of CR model in the form of list
  # x: the dose
  # a,b: parameter values of the CR model a=(a1,a2), b=(b1,b2)
  
  e1 <- exp(a[1]+b[1]*x)
  e2 <- exp(a[2]+b[2]*x)
  p1 <- 1/(1+e1) * 1/(1+e2) # ... p0 ... no reaction
  p2 <- e2 * p1 # ... pS ... success
  p3 <- e1 / (1+e1) # ... pT ... toxicity
  return(list(p1=p1,p2=p2,p3=p3))
}

resp2 <- function(x, a, b){
  #the response of CR model in the form of vector
  # x: the dose
  # a,b: parameter values of the CR model a=(a1,a2), b=(b1,b2)
  
  e1 <- exp(a[1]+b[1]*x)
  e2 <- exp(a[2]+b[2]*x)
  p1 <- 1/(1+e1) * 1/(1+e2) # ... p0 ... no reaction
  p2 <- e2 * p1 # ... pS ... success
  p3 <- e1 / (1+e1) # ... pT ... toxicity
  return(c(p1,p2,p3))
}

#------------------------------------------------------
#---------------- main function -----------------------
#------------------------------------------------------

CRLAS <- function(x, N, A, C, bk, a, b, type="exact", solver="MOSEK"){
  # function for computing constrained exact designs in CR model
  # x: set of doses (size-n vector)
  # N: size of the design (integer)
  # A: linear constraints (K*n matrix)
  # C: 'level of sparsity' constraints (K*n matrix)
  # bk: the right side of constraints Aw + C*sw <= bk (size-K vector)
  # a,b: parameter values of the CR model (size-2 vectors, a=(a1,a2), b=(b1,b2))
  
  
  lx <- length(x)
  
  Fx <- ContF(x[1], a, b)
  for (i in 2:lx) Fx <- rbind(Fx, ContF(x[i], a, b))
  
  # Re-ordering the regressors
  Fx.odd <- Fx[seq(1,2*lx, 2),]
  Fx.even <- Fx[seq(2, 2*lx, 2),]
  
  # Add zero regressors for y-s
  Fx <- rbind(Fx.odd, Fx.even, matrix(rep(0, 4*lx), nrow = lx, ncol = 4))
  
  
  #equality constraints for the multiresponse model
  # w(x_i,1) = w(x_i, 2) = ... = w(x_i, r) (5)
  A31 <- matrix(0, nrow=lx, ncol=3*lx)
  for (i in 1:lx){
    A31[i, i] <- 1
    A31[i, i + lx] <- -1
  } 
  b31 <- rep(0, lx)
  
  # sum w (xi, 1)   =  N (6)
  A32 <- c(rep(1, lx), rep(0, 2*lx))
  b32 <- N
  
  # w(zi) \in {0,1} (7)
  A11 <- cbind( matrix(rep(0, 2*lx^2), ncol = 2*lx, nrow = lx), diag(lx) )
  b11 <- rep(1, lx)
  
  # w(xi, 1) <= N*w(zi) (9)
  a12 <- c(1, -N)
  A12 <- matrix(0, nrow=lx, ncol=3*lx)
  for(i in 1:lx) A12[i, c(i, 2*lx+i)] <- a12
  b12 <- rep(0, lx)
  
  # w(zi) <= w(xi, 1) (9)
  a14 <- c(-1, 1)
  A14 <- matrix(0, nrow = lx, ncol = 3*lx)
  for(i in 1:lx) A14[i, c(i, 2*lx+i)] <- a14
  b14 <- rep(0, lx)
  
  # All w should be non-negative:
  A21 <- diag((3*lx))
  b21 <- rep(0, 3*lx)
  
  # the core: linear and sparsity constraints
  # K is the number of constraints
  K <- dim(A)[1] 
  A13 <- matrix(0, nrow=K, ncol=3*lx)
  A13[1:K, 1:lx] <- A
  A13[1:K, (2*lx+1):(3*lx)] <- C
  b13 <- bk
  
  #stitch it all together
  b1 <- c(b11, b12, b13, b14)
  A1 <- rbind(A11, A12, A13, A14)
  b2 <- b21
  A2 <- A21
  b3 <- c(b31, b32)
  A3 <- rbind(A31, A32)
  
  if (type=="approximate"){
    res <- MISOCP_m4(Fx, b1=b1, A1=A1, b2=b2, A2=A2, b3=b3, A3=A3, type=type, solver=solver)
  } else {
    res <- MISOCP_m4(Fx, b1=b1, A1=A1, b2=b2, A2=A2, b3=b3, A3=A3, type="exact", gap=0, solver=solver)
  }
    
  w.best <- rep(0, lx)
  for (i in 1:lx) w.best[i] <- res$w.best[i] #+ res$w.best[lx+i]
  supp <- (1:lx)[w.best > 1e-6]
  w.supp <- w.best[supp]
  xsupp <- x[supp]
  return(list(res=res, w.best=w.best, xsupp=xsupp, w.supp=w.supp)) 
}

#----------------------------------------------------------------
#--------------- Example in Section 5 ---------------------------
#----------------------------------------------------------------

# Initial parameter settings
a <- c(-9.5,-9.1)  # parameters of the CR model 
b <- c(0.12,0.33)  # parameters of the CR model
x <- seq(0, 100, length=101)  # vector of doses
n <- length(x)  # number of doses
N <- 100  # limit on the number of patients

type <- "exact"
crit <- "D"
ps <- resp(x, a, b)
pf <- ps$p1 + ps$p3  # probability of failure

ps2 <- resp2(x[1], a, b)
for(i in 2:n) ps2 <- rbind(ps2, resp2(x[i], a, b))

# Additional costs based on the output of the experiment

c0 <- 5 # no reaction
cs <- 0 # efficacy
ct <- 20 # toxicity

expected_costs <- t(as.matrix(ps2) %*% as.matrix(c(c0, cs, ct)))
expected_costs
overhead_costs <- 0.4*x # manufacturing costs 
overhead_costs

#### Size constraint only (N=100) ####

A <- matrix(0)
C <- matrix(0)
bk <- 0

res_0 <- CRLAS(x, N, A, C, bk, a, b)
Phi_0 <- res_0$res$Phi.best # for the efficiencies of constrained optimal designs
S_0 <- length(res_0$w.supp) # number of support points
Ef_0 <- t(as.matrix(pf)) %*% res_0$w.best # expected number of failures with no constraints
expected_costs_0 <- expected_costs %*% res_0$w.best # cost of the experiment 
overhead_costs_0 <- sum(overhead_costs %*% (res_0$w.best!=0)) 
total_costs_0 <- expected_costs_0 + overhead_costs_0

#### Expected number of failures, Sec. 5.1.1 ####

A1 <- t(as.matrix(pf))
C1 <- t(as.matrix(rep(0, n)))
bk1 <- 40

A <- A1
C <- C1
bk <- bk1

res_1 <- CRLAS(x, N, A, C, bk, a, b)

Phi_1 <- res_1$res$Phi.best
S_1 <- length(res_1$w.supp)
Ef_1 <- t(as.matrix(pf)) %*% res_1$w.best
expected_costs_1 <- expected_costs %*% res_1$w.best
overhead_costs_1 <- sum(overhead_costs %*% (res_1$w.best!=0)) 
total_costs_1 <- expected_costs_1 + overhead_costs_1

#### Constraint on the total costs of the experiment, Sec. 5.1.2. ####

A2 <- expected_costs
C2 <- overhead_costs
total_costs_limit <- 500
bk2 <- total_costs_limit

A <- rbind(A1, A2)
C <- rbind(C1, C2)
bk <- c(bk1, bk2)

res_2 <- CRLAS(x, N, A, C, bk, a, b)

Phi_2 <- res_2$res$Phi.best
S_2 <- length(res_2$w.supp)
Ef_2 <- t(as.matrix(pf)) %*% res_2$w.best
expected_costs_2 <- expected_costs %*% res_2$w.best
overhead_costs_2 <- sum(overhead_costs %*% (res_2$w.best!=0)) 
total_costs_2 <- expected_costs_2 + overhead_costs_2

#### Minimal support size, Sec. 5.1.3 ####

S <- 6
A3 <- t(rep(0,n))
C3 <- t(rep(-1, n))
bk3 <- -S

A <- rbind(A1, A2, A3)
C <- rbind(C1, C2, C3)
bk <- c(bk1, bk2, bk3)

res_3 <- CRLAS(x, N, A, C, bk, a, b)

Phi_3 <- res_3$res$Phi.best
S_3 <- length(res_3$w.supp)
Ef_3 <- t(as.matrix(pf)) %*% res_3$w.best
expected_costs_3 <- expected_costs %*% res_3$w.best
overhead_costs_3 <- sum(overhead_costs %*% (res_3$w.best!=0))
total_costs_3 <- expected_costs_3 + overhead_costs_3

#### Space filling, Sec. 5.1.4 ####

space <- 10 # the required spacing of the doses
spaces <- rep(1, space) # a sequence of ones
C4 <- matrix(0, nrow = n - space + 1, ncol = n)
for(i in 1:nrow(C4)){
  C4[i, i:(i+space-1)] <- spaces
}

A4 <- matrix(0, nrow = n - space + 1, ncol = n)

bk4 <- rep(1, n - space + 1)

A <- rbind(A1, A2, A3, A4)
C <- rbind(C1, C2, C3, C4)
bk <- c(bk1, bk2, bk3, bk4)

res_4 <- CRLAS(x, N, A, C, bk, a, b)

Phi_4 <- res_4$res$Phi.best
S_4 <- length(res_4$w.supp)
Ef_4 <- t(as.matrix(pf)) %*% res_4$w.best
expected_costs_4 <- expected_costs %*% res_4$w.best
overhead_costs_4 <- sum(overhead_costs %*% (res_4$w.best!=0)) 
total_costs_4 <- expected_costs_4 + overhead_costs_4

#### Minimal and maximal number of replications, Sec. 5.1.5. ####

L <- 10 # lower limit on the number of replications
U <- 25 # upper limit on the number of replications

A5 <- -diag(n)
C5 <- L*diag(n)
bk5 <- rep(0,n)

A6 <- diag(n)
C6 <- -U*diag(n)
bk6 <- rep(0,n)

A <- rbind(A1, A2, A3, A4, A5, A6)
C <- rbind(C1, C2, C3, C4, C5, C6)
bk <- c(bk1, bk2, bk3, bk4, bk5, bk6)

res_5 <- CRLAS(x, N, A, C, bk, a, b)

Phi_5 <- res_5$res$Phi.best
S_5 <- length(res_5$w.supp)
Ef_5 <- t(as.matrix(pf)) %*% res_5$w.best
expected_costs_5 <- expected_costs %*% res_5$w.best
overhead_costs_5 <- sum(overhead_costs %*% (res_5$w.best!=0)) 
total_costs_5 <- expected_costs_5 + overhead_costs_5

# Efficiencies:
c(Phi_1/Phi_0, Phi_2/Phi_0, Phi_3/Phi_0, Phi_4/Phi_0, Phi_5/Phi_0)


#---------- Generating Figure 1-----------------------------------#

library(ggplot2)
library(tidyr)
library(dplyr)

#0
l <- cbind(x,res_0$w.best)
l <- as.data.frame(l)
names(l) <- c("x", "w0*")
df <- pivot_longer(l, cols = -x, names_to = "design", values_to = "value")
df$x_jitter <- as.numeric(df$x) + (as.numeric(factor(df$design)) - 3) * 0.1
p0 <- ggplot(df, aes(x = x_jitter, xend = x_jitter, y = 0, yend = value, color = design)) +
  geom_segment(size = 2) +
  ylim(0, 45) +
  scale_x_continuous(breaks=l$x[l$x %% 5 == 0], labels = l$x[l$x %% 5 == 0]) +
  labs(x = "Dose", y = "Weight") +  theme_minimal()

p0 <- ggplot(df, aes(x = x_jitter, xend = x_jitter, y = 0, yend = value, color = design)) +
  geom_segment(size = 2) +
  ylim(0, 45) +
  scale_x_continuous(breaks = l$x[l$x %% 5 == 0], labels = l$x[l$x %% 5 == 0]) +
  labs(x = "Dose", y = "Weight") +
  ggtitle("Efficiency relative to w0*: 100%") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 1, size = 12, face = "plain"),  
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)              
  )

#1
l <- cbind(x,res_1$w.best)
l <- as.data.frame(l)
names(l) <- c("x", "w1*")
df <- pivot_longer(l, cols = -x, names_to = "design", values_to = "value")
df$x_jitter <- as.numeric(df$x) + (as.numeric(factor(df$design)) - 3) * 0.1
p1 <- ggplot(df, aes(x = x_jitter, xend = x_jitter, y = 0, yend = value, color = design)) +
  geom_segment(size = 2) +
  ylim(0, 45) +
  scale_x_continuous(breaks=l$x[l$x %% 5 == 0], labels = l$x[l$x %% 5 == 0]) +
  labs(x = "Dose", y = "Weight") +  theme_minimal()

p1 <- ggplot(df, aes(x = x_jitter, xend = x_jitter, y = 0, yend = value, color = design)) +
  geom_segment(size = 2) +
  ylim(0, 45) +
  scale_x_continuous(breaks = l$x[l$x %% 5 == 0], labels = l$x[l$x %% 5 == 0]) +
  labs(x = "Dose", y = "Weight") +
  ggtitle("Efficiency relative to w0*: 98%") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 1, size = 12, face = "plain"),  
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)              
  )
#2
l <- cbind(x,res_2$w.best)
l <- as.data.frame(l)
names(l) <- c("x", "w2*")
df <- pivot_longer(l, cols = -x, names_to = "design", values_to = "value")
df$x_jitter <- as.numeric(df$x) + (as.numeric(factor(df$design)) - 3) * 0.1
p2 <- ggplot(df, aes(x = x_jitter, xend = x_jitter, y = 0, yend = value, color = design)) +
  geom_segment(size = 2) +
  ylim(0, 45) +
  scale_x_continuous(breaks=l$x[l$x %% 5 == 0], labels = l$x[l$x %% 5 == 0]) +
  labs(x = "Dose", y = "Weight") +  theme_minimal()

p2 <- ggplot(df, aes(x = x_jitter, xend = x_jitter, y = 0, yend = value, color = design)) +
  geom_segment(size = 2) +
  ylim(0, 45) +
  scale_x_continuous(breaks = l$x[l$x %% 5 == 0], labels = l$x[l$x %% 5 == 0]) +
  labs(x = "Dose", y = "Weight") +
  ggtitle("Efficiency relative to w0*: 96%") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 1, size = 12, face = "plain"),  
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)              
  )

#3
l <- cbind(x,res_3$w.best)
l <- as.data.frame(l)
names(l) <- c("x", "w3*")
df <- pivot_longer(l, cols = -x, names_to = "design", values_to = "value")
df$x_jitter <- as.numeric(df$x) + (as.numeric(factor(df$design)) - 3) * 0.1
p3 <- ggplot(df, aes(x = x_jitter, xend = x_jitter, y = 0, yend = value, color = design)) +
  geom_segment(size = 2) +
  ylim(0, 45) +
  scale_x_continuous(breaks=l$x[l$x %% 5 == 0], labels = l$x[l$x %% 5 == 0]) +
  labs(x = "Dose", y = "Weight") +  theme_minimal()

p3 <- ggplot(df, aes(x = x_jitter, xend = x_jitter, y = 0, yend = value, color = design)) +
  geom_segment(size = 2) +
  ylim(0, 45) +
  scale_x_continuous(breaks = l$x[l$x %% 5 == 0], labels = l$x[l$x %% 5 == 0]) +
  labs(x = "Dose", y = "Weight") +
  ggtitle("Efficiency relative to w0*: 95%") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 1, size = 12, face = "plain"), 
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)               
  )

#4
l <- cbind(x,res_4$w.best)
l <- as.data.frame(l)
names(l) <- c("x", "w4*")
df <- pivot_longer(l, cols = -x, names_to = "design", values_to = "value")
df$x_jitter <- as.numeric(df$x) + (as.numeric(factor(df$design)) - 3) * 0.1
p4 <- ggplot(df, aes(x = x_jitter, xend = x_jitter, y = 0, yend = value, color = design)) +
  geom_segment(size = 2) +
  ylim(0, 45) +
  scale_x_continuous(breaks=l$x[l$x %% 5 == 0], labels = l$x[l$x %% 5 == 0]) +
  labs(x = "Dose", y = "Weight") +  theme_minimal()

p4 <- ggplot(df, aes(x = x_jitter, xend = x_jitter, y = 0, yend = value, color = design)) +
  geom_segment(size = 2) +
  ylim(0, 45) +
  scale_x_continuous(breaks = l$x[l$x %% 5 == 0], labels = l$x[l$x %% 5 == 0]) +
  labs(x = "Dose", y = "Weight") +
  ggtitle("Efficiency relative to w0*: 94%") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 1, size = 12, face = "plain"), 
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)               
  )

#5
l <- cbind(x,res_5$w.best)
l <- as.data.frame(l)
names(l) <- c("x", "w5*")
df <- pivot_longer(l, cols = -x, names_to = "design", values_to = "value")
df$x_jitter <- as.numeric(df$x) + (as.numeric(factor(df$design)) - 3) * 0.1
p5 <- ggplot(df, aes(x = x_jitter, xend = x_jitter, y = 0, yend = value, color = design)) +
  geom_segment(size = 2) +
  ylim(0, 45) +
  scale_x_continuous(breaks=l$x[l$x %% 5 == 0], labels = l$x[l$x %% 5 == 0]) +
  labs(x = "Dose", y = "Weight") +  theme_minimal()

p5 <- ggplot(df, aes(x = x_jitter, xend = x_jitter, y = 0, yend = value, color = design)) +
  geom_segment(size = 2) +
  ylim(0, 45) +
  scale_x_continuous(breaks = l$x[l$x %% 5 == 0], labels = l$x[l$x %% 5 == 0]) +
  labs(x = "Dose", y = "Weight") +
  ggtitle("Efficiency relative to w0*: 89%") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 1, size = 12, face = "plain"), 
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)      
  )

library(gridExtra)
grid.arrange(p0, p1, p2, p3, p4, p5, nrow=6)


