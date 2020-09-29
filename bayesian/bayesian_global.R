update_alpha.lasso <- function(datalist, alpha.old, sigma.square,theta, maxstep_inner, tolerance_inner) {
  ## initial
  alpha.inner.old = alpha.old
  k_inner = 1
  n = nrow(datalist[[1]][[1]])
  p = ncol(datalist[[1]][[1]])
  delta.est = list()
  while (k_inner < maxstep_inner) {
    # given alpha update delta

    delta.est = lapply(1:length(datalist), function(index){
      gamma = 2 * exp(-2 * datalist[[index]][[3]] %*% alpha.old)
      sd_y = sqrt(var(datalist[[index]][[2]]) * (n - 1)/n)
      C = sum(1/gamma)/p * sd_y * sigma.square[[index]]/n
      pen_vec = 1/gamma
      #pen_vec[which(!overlap_m[index, -index])] = Inf
      delta.est = coef(glmnet(datalist[[index]][[1]], datalist[[index]][[2]], alpha = 0, penalty.factor = pen_vec, lambda = C, standardize = FALSE, intercept = FALSE))[-1]
      return(delta.est)
    })

    # for(j in 1:length(datalist)){
    #   gamma = 2 * exp(-2 * datalist[[j]][[3]] %*% alpha.old)
    #   sd_y = sqrt(var(datalist[[j]][[2]]) * (n - 1)/n)
    #   C = sum(1/gamma)/p * sd_y * sigma.square/n
    #   pen_vec = 1/gamma
    #   pen_vec[which(!overlap_m[j, -j])] = Inf
    #   delta.est[[j]] = coef(glmnet(datalist[[j]][[1]], datalist[[j]][[2]], alpha = 0, penalty.factor = pen_vec, lambda = C, standardize = F, intercept = FALSE))[-1]
    # }
        ## given delta update alpha
    alpha.inner.new <- optim(alpha.old, likelihood.alpha.theta.lasso, likelihood.alpha.theta.gradient.lasso, datalist = datalist, theta = theta, delta = delta.est, method = "BFGS")$par
    if (sum(abs(alpha.inner.new - alpha.inner.old)) < tolerance_inner) {
      break
    }
    k_inner = k_inner + 1
    alpha.inner.old <- alpha.inner.new
  }
  return(list(alpha.est = alpha.inner.old, inner_iter = k_inner))
}

likelihood.alpha.theta.lasso <- function(datalist, alpha, theta, delta) {

  result = lapply(1:length(datalist), function(index){
    gamma = 2 * exp(-2 * datalist[[index]][[3]] %*% alpha)
    result = as.numeric(t(theta[[index]]) %*% gamma + delta[[index]]^2 %*% (1/gamma))
    return(result)
  })

#   result = rep(0, length(datalist))
#   for(i in 1: length(datalist)){
#     gamma = 2 * exp(-2 * datalist[[i]][[3]] %*% alpha)
#     result[i] = as.numeric(t(theta[[i]]) %*% gamma + delta[[i]]^2 %*% (1/gamma))
#   }
  return(sum(unlist(result)))
}

likelihood.alpha.theta.gradient.lasso <- function(datalist, alpha, theta, delta) {

  result = lapply(1:length(datalist), function(index){
      gamma = 2 * exp(-2 * datalist[[index]][[3]] %*% alpha)
      dev_gamma = (theta[[index]] - delta[[index]]^2/(gamma^2))
      result = crossprod(dev_gamma, as.vector(gamma) * datalist[[index]][[3]]) * (-2)
  })
  result = matrix(unlist(result), nrow = 5)
  # result = matrix(0, ncol = 4, nrow = length(datalist))
  # for(i in 1: length(datalist)){
  #   gamma = 2 * exp(-2 * datalist[[i]][[3]] %*% alpha)
  #   dev_gamma = (theta[[i]] - delta[[i]]^2/(gamma^2))
  #   result[i,] = crossprod(dev_gamma, as.vector(gamma) * datalist[[i]][[3]]) * (-2)
  # }
  #return(colSums(result))
  return(rowSums(result))
}

#############################################################
ptune <- function(docM, group, maxstep, tolerance, maxstep_inner, tolerance_inner){
  datalist = list()
  for(i in 1: ncol(docM)){
    y = docM[,i]
    logY = log(docM[,i]+1)
    X = docM[,-i]
    Z = lapply(group, function(m){return(m[-i,i])})
    Z = matrix(c(rep(1, ncol(X)),as.numeric(unlist(Z))), ncol = length(group)+1)
    datalist[[i]] = list(X, logY, Z)
  }
    q = length(group)
    alpha.old = as.matrix(rep(0, q+1))
    k = 1
    #sigma = estimateVariance(datalist[[1]][[1]], datalist[[1]][[2]])
    #sigma.square = rep(sigma, length(datalist))
    n = nrow(docM)
    sigma.square = lapply(datalist, function(datalist){return(mean(1/(datalist[[2]]+1)/n))})

    while (k < maxstep) {
      # Given alpha, update theta
      print(k)
      theta = lapply(1:length(datalist), function(index){
        gamma = 2 * exp(-2 * datalist[[index]][[3]] %*% alpha.old)
        Sigma_y = sigma.square[[index]] * diag(n) + (t(t(datalist[[index]][[1]]) * c(gamma))) %*% t(datalist[[index]][[1]])
        theta = colSums(datalist[[index]][[1]] * solve(Sigma_y, datalist[[index]][[1]]))
        return(theta)
      })

      # theta = list()
      # for(j in 1:length(datalist)){
      #   gamma = 2 * exp(-2 * datalist[[j]][[3]] %*% alpha.old)
      #   Sigma_y = diag(1/exp(datalist[[j]][[2]])) + (t(t(datalist[[j]][[1]]) * c(gamma))) %*% t(datalist[[j]][[1]])
      #   theta[[j]] = colSums(datalist[[j]][[1]] * solve(Sigma_y, datalist[[j]][[1]]))
      # }

      # Given theta, update alpha
      update.result <- update_alpha.lasso(datalist, alpha.old = alpha.old, sigma.square = sigma.square, theta = theta, maxstep_inner = maxstep_inner, tolerance_inner = tolerance_inner)
      alpha.new <- update.result$alpha.est

      # Check convergence
      if (sum(abs(alpha.new - alpha.old)) < tolerance) {
        print(k)
        break
      }
      alpha.old <- alpha.new

      k <- k + 1
    }

    return(list(alpha.est = alpha.old, sigma.square = sigma.square))
}
##########################################################
library(glmnet)
library(data.table)
library(igraph)
library(ggplot2)
#library(xtune)

load("data_community_large.Rdata")
load("penalties_large.Rdata")
people_list = data$ID

n = 2000
p = nrow(data)
set.seed(10)

generate_data = function(n, p, relation_list, persons){

  E = matrix(rpois(p*n, 0.004), nrow = n)
  Y = matrix(rpois(n*(p+p*(p-1)/2), 0.004), nrow = n)

  p1 = t(combn(1:p, 2))
  P = matrix(0, ncol = p, nrow = p*(p-1)/2)
  for(j in 1: (p*(p-1)/2)){
    P[j, p1[j, 1]] = 1
    P[j, p1[j, 2]] = 1
  }

  true_relation = matrix(rep(0, p*p), nrow = p)
  for(j in 1: nrow(relation_list)){
    ID1 = which(persons == relation_list$V1[j])
    ID2 = which(persons == relation_list$V2[j])
    true_relation[ID1, ID2] = 1
    true_relation[ID2, ID1] = 1
  }

  A = true_relation
  A[lower.tri(A, diag = TRUE)] <- NA
  A = as.vector(t(A))
  A = A[!is.na(A)]

  B = rbind(diag(p), P * t(rep(1, p) %*% t(A)))

  X = Y %*% B + E

  return(X)
}

doc_counts = generate_data((2*n), p, relationship, data$ID)
doc_counts = doc_counts[which(rowSums(doc_counts) != 0),]
doc_counts = doc_counts[1:n,]
overlap_m[which(overlap_m == FALSE)] = Inf
###################################################
start_time <- Sys.time()
alpha_p = ptune(doc_counts, group = list(family_m, group_A, group_B, group_C), maxstep = 100, tolerance = 0.01, maxstep_inner = 50, tolerance_inner = 0.1)
end_time <- Sys.time()
end_time - start_time

save(alpha_p, file = "alpha10.Rdata")