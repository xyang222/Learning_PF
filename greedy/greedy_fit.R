library(glmnet)
library(data.table)
library(igraph)
library(ggplot2)

load("data_community_large.Rdata")
load("penalties_large.Rdata")
people_list = data$ID

################################################simu data
#600 and 2000
n = 2000
p = nrow(data)

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

################################################
overlap_m[which(overlap_m == FALSE)] = Inf
overlap_m[which(overlap_m == TRUE)] = 0

alpha_family = seq(-1, -0.5, 0.05)
alpha_group_A = seq(-0.5, 0, 0.05)
alpha_group_B = seq(-0.5, 0, 0.05)
alpha_group_C = seq(-0.5, 0, 0.05)

#alpha = expand.grid(alpha_family, alpha_group)
#alpha = expand.grid(alpha_family, alpha_group_A, alpha_group_B, alpha_group_C)

#group = (group_A|group_B)

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

fit_cv <- function(person_counts, other_counts, penalty_list) {
  lambdas = seq(0.005, 0.1, 0.01)
  fit.cv = cv.glmnet(other_counts, person_counts, family='poisson', thresh = 1e-06, lambda = lambdas, maxit = 2000, type.measure="mse", penalty.factor = penalty_list, lower.limits = 0)
  return(fit.cv)
}

local_poisson = function(doc_counts, penalty, persons){

  results = lapply(1:length(persons), function(i) {
    y = doc_counts[,i]
    X = doc_counts[,-i]
    penalty_list = penalty[-i,i]
    result = fit_cv(y, X, penalty_list)
    return(result)
  })
  nodes = as.data.frame(persons)

  mse = sapply(1:length(results), function(i){
    result = results[i][[1]]
    se.min <- result$cvm[result$lambda == result$lambda.min]
  })

  return(sum(mse))

}

fit_model = function(alpha, doc_counts, persons){

  print(alpha)

  final_penalty = exp(overlap_m*1 + alpha[1] * (0+family_m) + alpha[2] * (0+group_A) + alpha[3] * (0+group_B) + alpha[4] * (0+group_C))

  lp_result = tryCatch(local_poisson(doc_counts, final_penalty, people_list), error = function(e){"Error"})

  if(lp_result[1] == "Error" ) return(c(1000))

  print(lp_result)
  return(lp_result)
}

model_all = function(index, alpha_family, alpha_group_A, alpha_group_B, alpha_group_C){

  set.seed(index)
  doc_counts = generate_data((2*n), p, relationship, data$ID)
  doc_counts = doc_counts[which(rowSums(doc_counts) != 0),]
  doc_counts = doc_counts[1:n,]

  initial = c(0,0,0,0)

  alpha1 = expand.grid(alpha_family, initial[2], initial[3], initial[4])

  result = apply(alpha1, 1, fit_model, doc_counts = doc_counts, persons = people_list)
  result = unlist(result)
  maxIndex1 = which(result == min(result))[1]

  alpha2 = expand.grid(alpha_family[maxIndex1], alpha_group_A,initial[3], initial[4])

  result = apply(alpha2, 1, fit_model, doc_counts = doc_counts, persons = people_list)
  result = unlist(result)
  maxIndex2 = which(result == min(result))[1]

  alpha3 = expand.grid(alpha_family[maxIndex1], alpha_group_A[maxIndex2], alpha_group_B, initial[4])

  result = apply(alpha3, 1, fit_model, doc_counts = doc_counts, persons = people_list)
  result = unlist(result)
  maxIndex3 = which(result == min(result))[1]

  alpha4 = expand.grid(alpha_family[maxIndex1], alpha_group_A[maxIndex2], alpha_group_B[maxIndex3], alpha_group_C)

  result = apply(alpha4, 1, fit_model, doc_counts = doc_counts, persons = people_list)
  result = unlist(result)
  maxIndex4 = which(result == min(result))[1]

  print(c(alpha_family[maxIndex1], alpha_group_A[maxIndex2], alpha_group_B[maxIndex3], alpha_group_C[maxIndex4]))

  return(c(alpha_family[maxIndex1], alpha_group_A[maxIndex2], alpha_group_B[maxIndex3], alpha_group_C[maxIndex4]))
}

start_time <- Sys.time()
result_all = lapply(1:5, model_all, alpha_family = alpha_family, alpha_group_A = alpha_group_A, alpha_group_B = alpha_group_B, alpha_group_C = alpha_group_C)
end_time <- Sys.time()
end_time - start_time

save(result_all, file = "alpha_greedy11.Rdata")

start_time <- Sys.time()
test5 = lapply(5, model_all, alpha_family = alpha_family, alpha_group_A = alpha_group_A, alpha_group_B = alpha_group_B, alpha_group_C = alpha_group_C)
end_time <- Sys.time()
end_time - start_time

#######################################################################
load("alpha_greedy1.Rdata")

alpha_list = matrix(unlist(result_all), byrow = TRUE, nrow = 5)

fit_cv <- function(person_counts, other_counts, penalty_list) {
  fit.cv = cv.glmnet(other_counts, person_counts, nlambda = 10, family = 'poisson', thresh = 1e-06, maxit = 2000, type.measure="mse", penalty.factor = penalty_list, lower.limits = 0)
  return(fit.cv)
}

local_poisson = function(doc_counts, penalty, persons){

  results = lapply(1:length(persons), function(i) {
    #print(i)
    y = doc_counts[,i]
    X = doc_counts[,-i]
    penalty_list = penalty[-i,i]
    result = fit_cv(y, X, penalty_list)
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

generate_result = function(index){
  set.seed(index)
  doc_counts = generate_data((2*n), p, relationship, data$ID)
  doc_counts = doc_counts[which(rowSums(doc_counts) != 0),]
  doc_counts = doc_counts[1:n,]

  alpha = alpha_list[index,]

  result = fit_model(alpha, doc_counts, people_list)
  return(result)
}

start_time <- Sys.time()
result1 = sapply(1:5, generate_result)
end_time <- Sys.time()
end_time - start_time

save(result1, file = "result_greedy11.Rdata")
