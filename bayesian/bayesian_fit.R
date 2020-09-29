library(glmnet)
library(data.table)
library(igraph)
library(ggplot2)
library(plyr)

setwd("C:/Users/Xiaoyi/Dropbox/PhD/Bacon_afterProposal/cross_test")
load("data_community_large.Rdata")
load("penalties_large.Rdata")
load("alpha2.Rdata")
people_list = data$ID

################################################simu data
#600 and 2000
n = 2000
p = nrow(data)
set.seed(2)

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

################################################
overlap_m[which(overlap_m == FALSE)] = Inf

################################################ fit model functions
evaluation_simu = function(links, relationList){

  if(nrow(links) == 0 & nrow(relationList) == 0)
    return(c(1,1))
  if(nrow(links) == 0 & nrow(relationList) != 0)
    return(c(1,0))
  if(nrow(links) != 0 & nrow(relationList) == 0)
    return(c(0,1))
  if(ncol(links) == 1){
    links = as.data.frame(t(as.matrix(links)))
  }

  links = as.data.frame(links)
  colnames(links) = c("V1", "V2")
  links$V1 = as.integer(as.character(links$V1))
  links$V2 = as.integer(as.character(links$V2))

  testp = fintersect(setDT(relationList), setDT(links[,1:2]))
  colnames(links) = c("V2", "V1")
  testn = fintersect(setDT(relationList[,2:1]), setDT(links[,1:2]))

  precision = (nrow(testp) +nrow(testn))/nrow(links)
  recall = (nrow(testp) +nrow(testn))/nrow(relationList)

  return(c(precision, recall))
}

fit_cv <- function(person_counts, other_counts, penalty_list, lower) {
  fit.cv = cv.glmnet(other_counts, person_counts, family='poisson',thresh = 1e-06, nlambda = 10, maxit = 2000, type.measure="mse", penalty.factor = penalty_list, lower.limits = lower)
  return(fit.cv)
}

local_poisson = function(doc_counts, alpha_est, sigma, persons, lower){

  results = lapply(1:length(persons), function(i) {
    print(i)
    y = doc_counts[,i]
    X = doc_counts[,-i]
    Z = as.matrix(data.frame(rep(1, ncol(X)), family_m[-i,i], group_A[-i,i], group_B[-i,i], group_C[-i,i]))
    pen_vec_1 = exp(Z %*% alpha_est)*sigma[i]/n
    pen_vec_1 = pen_vec_1 * overlap_m[i, -i]
    result = fit_cv(y, X, pen_vec_1, lower)
    return(result)
  })
  nodes = as.data.frame(persons)

  links = lapply(1:length(results), function(i) {
      result = results[i][[1]] # Unlist it
      node = persons[i]
      coefs = coef(result, s = "lambda.min")
      coefs = coefs[-1] # Remove intercept
      relationships_idx = which(coefs > 0)
      relationships = persons[-i][relationships_idx] # Need to first remove the 'target' from the list so that the other idxs line up
      relationship_weights = coefs[relationships_idx]
      if(length(relationships) > 0){
        result.df = data.frame(a=node, b=relationships, w=relationship_weights)
      }else{
        result.df = data.frame()
      }
      return(result.df)
    })
  links.df = rbindlist(links)

    if(nrow(links.df) > 0){
      mydf = links.df[,1:2]
      links.final = links.df[!duplicated(apply(mydf,1,function(x) paste(sort(x),collapse=''))),]
    }else{
      links.final = data.frame()
    }

  result_l = evaluation_simu(links.final, relationship)

  return(c(result_l, nrow(links.final)))

}

model_fit = function(alpha){

  print(alpha)

  #lp_result = local_poisson(doc_counts, final_penalty, people_list)
  lp_result = tryCatch(local_poisson(doc_counts, alpha, people_list, lower = 0), error = function(e){"Error"})

  if(lp_result[1] == "Error" ) return(c(0,0,0))

  lp_maxIndex = which(rowSums(lp_result[,1:2]) == max(rowSums(lp_result[,1:2])))
  lp_final = lp_result[lp_maxIndex[1],]
  print(lp_final)

  return(lp_final)
}

local_poisson_g = function(doc_counts, penalty, persons, lower){

  results = lapply(1:length(persons), function(i) {
    #print(i)
    y = doc_counts[,i]
    X = doc_counts[,-i]
    penalty_list = penalty[-i,i]
    result = fit_cv(y, X, penalty_list, lower)
    return(result)
  })
  nodes = as.data.frame(persons)

  #paraList = c(seq(0.005, 0.009, 0.001) ,seq(0.01, 0.05, 0.01))
  paraList = seq(0.005, 0.05, 0.005)
  #paraList = seq(0.001, 0.05, 0.001)
  l_result = matrix(0, nrow = length(paraList), ncol = 3)
  for(k in 1:length(paraList)){
    links = lapply(1:length(results), function(i) {
      result = results[i][[1]] # Unlist it
      node = persons[i]
      coefs = coef(result, s = paraList[k])
      coefs = coefs[-1] # Remove intercept
      relationships_idx = which(coefs > 0)
      relationships = persons[-i][relationships_idx] # Need to first remove the 'target' from the list so that the other idxs line up
      relationship_weights = coefs[relationships_idx]
      if(length(relationships) > 0){
        result.df = data.frame(a=node, b=relationships, w=relationship_weights)
      }else{
        result.df = data.frame()
      }
      return(result.df)
    })
    links.df = rbindlist(links)

    if(nrow(links.df) > 0){
      mydf = links.df[,1:2]
      links.final = links.df[!duplicated(apply(mydf,1,function(x) paste(sort(x),collapse=''))),]
    }else{
      links.final = data.frame()
    }

    result_l = evaluation_simu(links.final, relationship)
    l_result[k,1:2] = result_l
    l_result[k, 3] = nrow(links.final)
  }
  return(l_result)

}

fit_model_g = function(alpha, doc_counts, persons){

  print(alpha)

  final_penalty = overlap_m*1 + alpha[1] * (1-family_m) + alpha[2] * (1-group_A) + alpha[3] * (1-group_B) + alpha[4] * (1-group_C)

  lp_result = tryCatch(local_poisson(doc_counts, final_penalty, people_list, lower = 0), error = function(e){"Error"})

  if(lp_result[1] == "Error" ) return(c(0,0,0))

  lp_maxIndex = which(rowSums(lp_result[,1:2]) == max(rowSums(lp_result[,1:2])))
  lp_final = lp_result[lp_maxIndex[1],]
  print(lp_final)
  return(lp_final)
}

generate_result = function(index){
  set.seed(index)
  doc_counts = generate_data((2*n), p, relationship, data$ID)
  doc_counts = doc_counts[which(rowSums(doc_counts) != 0),]
  doc_counts = doc_counts[1:n,]

  filename = paste("alpha", index, ".Rdata", sep = "")
  load(filename)

  result_b = local_poisson(doc_counts, alpha_p$alpha.est, unlist(alpha_p$sigma.square), people_list, 0)
  result_s = local_poisson(doc_counts, c(0,0,0,0,0), unlist(alpha_p$sigma.square), people_list, 0)

  return(c(result_b, result_s))

}

###############################################################

model_fit(alpha_p)

result = apply(alpha, 1, model_fit)
#save(result, file = "result_community_large.Rdata")

####################################################

result = sapply(1:10, generate_result)
save(result, file = "result_bayesian.Rdata")

load("result_bayesian.Rdata")

result_data = data.frame(data = c(result[1,], result[4,], result[2,], result[5,]), group = c(rep("precision", 20), rep("recall", 20)), method = c(rep("Bayesian", 10), rep("SDFB", 10), rep("Bayesian", 10), rep("SDFB", 10)))

load("result_greedy11.Rdata")
greedy1 = data.frame(data = c(result1[1,], result1[2,]), group = c(rep("precision", 5), rep("recall", 5)), method = rep("Greedy", 10))
load("result_greedy21.Rdata")
greedy2 = data.frame(data = c(result2[1,], result2[2,]), group = c(rep("precision", 5), rep("recall", 5)), method = rep("Greedy", 10))

result_data = rbind(result_data, greedy1, greedy2)

average = data.frame(data = c((c(result1[1,], result2[1,])+c(result1[2,], result2[2,]))/2, (result[1,]+result[2,])/2, (result[4,]+result[5,])/2), group = c(rep("average", 30)), method = c(rep("Greedy", 10), rep("Bayesian", 10), rep("SDFB", 10)))

result_data = rbind(result_data, average)

g = ggplot(result_data, aes(y = data, x = factor(method, levels = c("SDFB", "Greedy", "Bayesian"))))+ geom_boxplot(data = result_data,  aes(group = method)) + facet_wrap(~factor(group, levels = c("precision", "recall", "average")),scales="free_x") + theme(axis.text=element_text(size=12), axis.title = element_text(size = 14), legend.text=element_text(size = 12), plot.title = element_text(size = 14)) + labs(x = "", y = "value")
plot(g)

##################################################################################

load("alpha_greedy11.Rdata")
alpha1 = matrix(unlist(result_all), byrow = TRUE, nrow = 5)
load("alpha_greedy21.Rdata")
alpha2 = matrix(unlist(result_all), byrow = TRUE, nrow = 5)

alpha_data = rbind(alpha1, alpha2)
x = alpha_data[,1]
y = 1:10
par(mfcol = c(2, 1))
plot(x,y,xlim = c(0, 1),ylim = c(1, 10), col = "white", xlab = "alpha", ylab = "runs", main = "Absolute value of alpha for Greedy method")
color = c("darkgreen", "black", "blue", "red")
pointshape = c(15,16,18,17)
for(j in 1:4){
  x = abs(alpha_data[,j])
  y = 1:10
  points(x,y, col = color[j], cex = 0.8, pch = pointshape[j])
}
#lines(x = seq(-0.011, 1.3, 0.001), y=rep(15.5, 1312), lty = 2)
legend("bottomright", legend = c("alpha_family", "alpha_groupA", "alpha_groupB", "alpha_groupC"), col = c("darkgreen", "black", "blue", "red"), pch = c(15,16,18,17), bty = "n", pt.cex = 0.8, cex = 0.8, text.col = "black", horiz = F , inset = c(0, 0))

alpha_data = matrix(0, nrow = 10, ncol = 4)
for(i in 1:10){
  filename = paste("alpha", i, ".Rdata", sep = "")
  load(filename)
  alpha_data[i,] = alpha_p$alpha.est[2:5]
}

x = alpha_data[,1]
y = 1:10
plot(x,y,xlim = c(0, 1),ylim = c(1, 10), col = "white", xlab = "alpha", ylab = "runs", main = "Absolute value of alpha for Bayesian method")
color = c("darkgreen", "black", "blue", "red")
pointshape = c(15,16,18,17)
for(j in 1:4){
  x = abs(alpha_data[,j])
  y = 1:10
  points(x,y, col = color[j], cex = 0.8, pch = pointshape[j])
}
#lines(x = seq(-0.011, 1.3, 0.001), y=rep(15.5, 1312), lty = 2)
legend("bottomright", legend = c("alpha_family", "alpha_groupA", "alpha_groupB", "alpha_groupC"), col = c("darkgreen", "black", "blue", "red"), pch = c(15,16,18,17), bty = "n", pt.cex = 0.8, cex = 0.8, text.col = "black", horiz = F , inset = c(0, 0))

