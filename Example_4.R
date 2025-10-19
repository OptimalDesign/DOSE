library(OptimalDesign)

#---------------- helper functions ---------------------#

ContF <- function(x, a, b){
  # elementary information matrices for continuation ratio (CR) model in dose x
  # parameter values a=c(a1,a2), b=c(b1,b2)
  
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
  #the response of CR model in the form of list
  e1 <- exp(a[1]+b[1]*x)
  e2 <- exp(a[2]+b[2]*x)
  p1 <- 1/(1+e1) * 1/(1+e2) # ... p0 ... no reaction
  p2 <- e2 * p1 # ... pS ... success
  p3 <- e1 / (1+e1) # ... pT ... toxicity
  return(list(p1=p1,p2=p2,p3=p3))
}

resp2 <- function(x, a, b){
  #the response of CR model in the form of vector
  e1 <- exp(a[1]+b[1]*x)
  e2 <- exp(a[2]+b[2]*x)
  p1 <- 1/(1+e1) * 1/(1+e2) # ... p0 ... no reaction
  p2 <- e2 * p1 # ... pS ... success
  p3 <- e1 / (1+e1) # ... pT ... toxicity
  return(c(p1,p2,p3))
}

CRLAS <- function(x, N, A, C, bk, a, b, crit="D", type="exact", w0=NULL){
  # function for computing constrained exact designs in CR model
  # x: set of doses
  # N: size of the design
  # A: linear constraints ... K*n matrix
  # C: 'level of sparsity' constraints ... K*n matrix
  # bk: the ride side of constraints Aw + C*sw <= bk
  # a,b: parameter values
  # crit: optimality criterion
  # type: "exact" or "approximate"
  # w0: design to be augmented
  
  lx <- length(x)
  
  Fx <- ContF(x[1], a, b)
  for (i in 2:lx) Fx <- rbind(Fx, ContF(x[i], a, b))
  
  # Re-ordering the regressors
  Fx.odd <- Fx[seq(1,2*lx, 2),]
  Fx.even <- Fx[seq(2, 2*lx, 2),]
  
  # Add zero regressors for y-s
  Fx <- rbind(Fx.odd, Fx.even, matrix(rep(0, 4*lx), nrow = lx, ncol = 4))
  
  #for (i in 1:lx) Fx[(4*i-1):(4*i),] <- Fx[(4*i-3):(4*i-2),]
  
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
    res <- od_MISOCP(Fx, b1=b1, A1=A1, b2=b2, A2=A2, b3=b3, A3=A3, type=type, crit=crit, w0=NULL)
  } else {
    res <- od_MISOCP(Fx, b1=b1, A1=A1, b2=b2, A2=A2, b3=b3, A3=A3, type="exact", crit=crit, w0=NULL, gap=0)
  }
  
  w.best <- rep(0, lx)
  for (i in 1:lx) w.best[i] <- res$w.best[i] #+ res$w.best[lx+i]
  supp <- (1:lx)[w.best > 1e-6]
  w.supp <- w.best[supp]
  xsupp <- x[supp]
  return(list(res=res, w.best=w.best, xsupp=xsupp, w.supp=w.supp)) 
}

#--------------- Example in Chapter 4--------------------------#

# Initial parameter settings
a <- c(-9.5,-9.1)  #parameters of the CR model 
b <- c(0.12,0.33)  #parameters of the CR model
x <- seq(0, 100, length=101)  #doses
n <- length(x)  #number of doses
N <- 100  #limit on the number of patients

type <- "exact"
crit <- "D"
ps <- resp(x, a, b)
pf <- ps$p1 + ps$p3 # probability of failure

ps2 <- resp2(x[1], a, b)
for(i in 2:n){
  ps2 <- rbind(ps2, resp2(x[i], a, b))
}

# Additional costs based on the output of the experiment

c0 <- 5 # no reaction
cs <- 0 # efficacy
ct <- 20 # toxicity

expected_costs <- t(as.matrix(ps2) %*% as.matrix(c(c0, cs, ct)))
expected_costs
overhead_costs <- 0.4*x # manufacturing costs 
overhead_costs

#### No constraints, D-optimality ####

A <- matrix(0)
C <- matrix(0)
bk <- 0

res_0 <- CRLAS(x, N, A, C, bk, a, b, crit="D", type="exact", w0=NULL)
Phi_0 <- res_0$res$Phi.best # for the efficiencies of constrained optimal designs
S_0 <- length(res_0$w.supp) # number of support points
Ef_0 <- t(as.matrix(pf)) %*% res_0$w.best # expected number of failures with no constraints
expected_costs_0 <- expected_costs %*% res_0$w.best # cost of the experiment 
overhead_costs_0 <- sum(overhead_costs %*% (res_0$w.best!=0)) 
total_costs_0 <- expected_costs_0 + overhead_costs_0

#### Expected number of failures ####

A1 <- t(as.matrix(pf))
C1 <- t(as.matrix(rep(0, n)))
bk1 <- 40

A <- A1
C <- C1
bk <- bk1

res_1 <- CRLAS(x, N, A, C, bk, a, b, crit="D", type="exact", w0=NULL)

Phi_1 <- res_1$res$Phi.best
S_1 <- length(res_1$w.supp)
Ef_1 <- t(as.matrix(pf)) %*% res_1$w.best
expected_costs_1 <- expected_costs %*% res_1$w.best
overhead_costs_1 <- sum(overhead_costs %*% (res_1$w.best!=0)) 
total_costs_1 <- expected_costs_1 + overhead_costs_1

#### Constraint on the total costs of the experiment

A2 <- expected_costs
C2 <- overhead_costs
total_costs_limit <- 500
bk2 <- total_costs_limit

A <- rbind(A1, A2)
C <- rbind(C1, C2)
bk <- c(bk1, bk2)

res_2 <- CRLAS(x, N, A, C, bk, a, b, crit="D", type="exact", w0=NULL)

Phi_2 <- res_2$res$Phi.best
S_2 <- length(res_2$w.supp)
Ef_2 <- t(as.matrix(pf)) %*% res_2$w.best
expected_costs_2 <- expected_costs %*% res_2$w.best
overhead_costs_2 <- sum(overhead_costs %*% (res_2$w.best!=0)) 
total_costs_2 <- expected_costs_2 + overhead_costs_2

#### Minimal support size ####

S <- 6
A3 <- t(rep(0,n))
C3 <- t(rep(-1, n))
bk3 <- -S

A <- rbind(A1, A2, A3)
C <- rbind(C1, C2, C3)
bk <- c(bk1, bk2, bk3)

res_3 <- CRLAS(x, N, A, C, bk, a, b, crit="D", type="exact", w0=NULL)

Phi_3 <- res_3$res$Phi.best
S_3 <- length(res_3$w.supp)
Ef_3 <- t(as.matrix(pf)) %*% res_3$w.best
expected_costs_3 <- expected_costs %*% res_3$w.best
overhead_costs_3 <- sum(overhead_costs %*% (res_3$w.best!=0))
total_costs_3 <- expected_costs_3 + overhead_costs_3

#### Space filling ####

space <- 10

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

res_4 <- CRLAS(x, N, A, C, bk, a, b, crit="D", type="exact", w0=NULL)

Phi_4 <- res_4$res$Phi.best
S_4 <- length(res_4$w.supp)
Ef_4 <- t(as.matrix(pf)) %*% res_4$w.best
expected_costs_4 <- expected_costs %*% res_4$w.best
overhead_costs_4 <- sum(overhead_costs %*% (res_4$w.best!=0)) 
total_costs_4 <- expected_costs_4 + overhead_costs_4

#### Minimal and maximal number of replications ####

L <- 10
U <- 25

A5 <- -diag(n)
C5 <- L*diag(n)
bk5 <- rep(0,n)

A6 <- diag(n)
C6 <- -U*diag(n)
bk6 <- rep(0,n)

A <- rbind(A1, A2, A3, A4, A5, A6)
C <- rbind(C1, C2, C3, C4, C5, C6)
bk <- c(bk1, bk2, bk3, bk4, bk5, bk6)

res_5 <- CRLAS(x, N, A, C, bk, a, b, crit="D", type="exact", w0=NULL)

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




