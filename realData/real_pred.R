c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A","#FF7F00", "black", "gold1","skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3", "darkorange4", "brown")

############################################################
fit_cv <- function(person_counts, other_counts, penalty_list) {
  fit.cv = cv.glmnet(other_counts, person_counts, family = 'poisson', thresh = 1e-06, maxit = 2000, type.measure="mse", penalty.factor = penalty_list, lower = 0)
  return(fit.cv)
}

get_links_bayesian = function(doc_counts, alpha_est, sigma, persons){

  results = lapply(1:length(persons), function(i) {
    #print(i)
    y = doc_counts[,i]
    X = doc_counts[,-i]
    Z = as.matrix(data.frame(rep(1, ncol(X)), family_m[-i,i], group_poet[-i,i], group_church[-i,i], group_royal[-i,i]))
    alpha_est = matrix(alpha_est, nrow = 5)
    pen_vec_1 = exp(Z %*% alpha_est)*sigma[i]/n
    pen_vec_1 = pen_vec_1 * overlap_m[i, -i]
    result = tryCatch(fit_cv(y, X, pen_vec_1), error = function(e){"Error"})
    if(result[1] == "Error" ) return(NULL)
    #result = fit_cv(y, X, pen_vec_1)
    return(result)
  })
  nodes = as.data.frame(persons)

  links = lapply(1:length(results), function(i) {
    result = results[i][[1]] # Unlist it
    if(is.null(result)) return(data.frame())
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

  return(links.final)

}

get_links_SDFB = function(doc_counts, persons){

  results = lapply(1:length(persons), function(i) {
    #print(i)
    y = doc_counts[,i]
    X = doc_counts[,-i]
    #Z = as.matrix(data.frame(rep(1, ncol(X)), family_m[-i,i], group_poet[-i,i], group_church[-i,i], group_royal[-i,i]))
    #alpha_est = matrix(alpha_est, nrow = 5)
    #pen_vec_1 = exp(Z %*% alpha_est)*sigma[i]/n
    pen_vec_1 = overlap_m[i, -i]
    result = tryCatch(fit_cv(y, X, pen_vec_1), error = function(e){"Error"})
    if(result[1] == "Error" ) return(NULL)
    #result = fit_cv(y, X, pen_vec_1)
    return(result)
  })
  nodes = as.data.frame(persons)

  links = lapply(1:length(results), function(i) {
    result = results[i][[1]] # Unlist it
    if(is.null(result)) return(data.frame())
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

  return(links.final)

}

bayesian = get_links_bayesian(doc_counts, alpha_p$alpha.est, unlist(alpha_p$sigma.square), people_list)
#SDFB = get_links_bayesian(doc_counts, c(0,0,0,0,0), unlist(alpha_p$sigma.square), people_list)
SDFB = get_links_SDFB(doc_counts, people_list)

save(bayesian, SDFB, file = "reallink1.Rdata")

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

truth = relationship[,2:3]
colnames(truth) = c("V1", "V2")
b_evalu = evaluation_simu(bayesian, truth)
s_evalu = evaluation_simu(SDFB, truth)

switchorder = function(data){
  for(i in 1:nrow(data)){
    if(data[i, 1] > data[i, 2]){
      temp = data[i, 2]
      data[i, 2] = data[i, 1]
      data[i, 1] = temp
    }
  }
  return(data)
}
b_switch = switchorder(bayesian[,1:2])
s_switch = switchorder(SDFB[,1:2])
t_switch = switchorder(truth[,1:2])
colnames(t_switch) = c("a", "b")
testdiff_b = fsetdiff(setDT(b_switch), setDT(s_switch))
testdiff_s = fsetdiff(setDT(s_switch), setDT(b_switch))

#######################################################

subdataID = unique(c(testdiff_b$a, testdiff_b$b, testdiff_s$a, testdiff_s$b))
subdata = people[which(people$SDFB_ID %in% subdataID),]

g_b = graph_from_data_frame(as.data.frame(testdiff_b), directed = FALSE, vertices = subdataID)
g_s = graph_from_data_frame(as.data.frame(testdiff_s), directed = FALSE, vertices = subdataID)

subrelation = relationship[which((relationship$person1_index %in% subdataID) & (relationship$person2_index %in% subdataID)),]
g = graph_from_data_frame(subrelation[,2:3], directed = FALSE, vertices = subdataID)
layout = layout_with_fr(g)

par(mfrow=c(1,2))
plot(g_b, vertex.label=NA, layout = layout, vertex.size=7)
plot(g_s, vertex.label=NA, layout = layout, vertex.size=7)

###########################################################
sort(table(c(testdiff_b$a, testdiff_b$b))) #10001185 10013000
sort(table(c(testdiff_s$a, testdiff_s$b)))
#10010449 10001346 10010925

#Robert Wingfield 10013000
link_b = c(testdiff_b$a[which(testdiff_b$b == 10013000)], testdiff_b$b[which(testdiff_b$a == 10013000)])
link_s = c(testdiff_s$a[which(testdiff_s$b == 10013000)], testdiff_s$b[which(testdiff_s$a == 10013000)])
people[which(people$SDFB_ID %in% link_b),]

#Elizabeth Blount 10001185
link_b = c(testdiff_b$a[which(testdiff_b$b == 10001185)], testdiff_b$b[which(testdiff_b$a == 10001185)])
link_s = c(testdiff_s$a[which(testdiff_s$b == 10001185)], testdiff_s$b[which(testdiff_s$a == 10001185)])
people[which(people$SDFB_ID %in% link_b),]

# Katherine Seymour 10010925
link_b = c(testdiff_b$a[which(testdiff_b$b == 10010925)], testdiff_b$b[which(testdiff_b$a == 10010925)])
link_s = c(testdiff_s$a[which(testdiff_s$b == 10010925)], testdiff_s$b[which(testdiff_s$a == 10010925)])
people[which(people$SDFB_ID %in% link_s),]

# Nicholas Bourbon 10001346
link_b = c(testdiff_b$a[which(testdiff_b$b == 10001346)], testdiff_b$b[which(testdiff_b$a == 10001346)])
link_s = c(testdiff_s$a[which(testdiff_s$b == 10001346)], testdiff_s$b[which(testdiff_s$a == 10001346)])
people[which(people$SDFB_ID %in% link_s),]

# Margaret Roper 10010449
$b[which(testdiff_b$a == 10010449)])
link_s = c(testdiff_s$a[which(testdiff_s$b == 10010449)], testdiff_s$b[which(testdiff_s$a == 10010449)])
people[which(people$SDFB_ID %in% link_s),]

################################################################

b_name1 = people$search_name1[match(testdiff_b$a, people$SDFB_ID)]
b_name2 = people$search_name1[match(testdiff_b$b, people$SDFB_ID)]
b_name = data.frame(b_name1, b_name2)
#write.csv(b_name, file = "bname.csv")

s_name1 = people$search_name1[match(testdiff_s$a, people$SDFB_ID)]
s_name2 = people$search_name1[match(testdiff_s$b, people$SDFB_ID)]
s_name = data.frame(s_name1, s_name2)
#write.csv(s_name, file = "sname.csv")

b_result = read.csv("bname.csv")
s_result = read.csv("sname.csv")
table(b_result$Relation)
#23/36  64.89%
table(s_result$Relation)
#28/59  47.46%
################################################################
subrelation3 = bayesian[which(bayesian$a < 1100 & bayesian$b < 1100),1:2]
colnames(subrelation3) = c("a", "b")

mydf = rbind(subrelation, subrelation3[,1:2])
union3 = mydf[!duplicated(apply(mydf,1,function(x) paste(sort(x),collapse=''))),]

color3 = rep("red", nrow(union3))
for(i in 1:nrow(union3)){
  index = which((subrelation3[,1] == as.numeric(union3[i,1]) & subrelation3[,2] == as.numeric(union3[i,2])) | (subrelation3[,1] == as.numeric(union3[i,2]) & subrelation3[,2] == as.numeric(union3[i,1])))
  if(length(index) > 0){
    color1[i] = "grey"
  }
  index = which((subrelation[,1] == as.numeric(union3[i,1]) & subrelation[,2] == as.numeric(union3[i,2])) | (subrelation[,1] == as.numeric(union3[i,2]) & subrelation[,2] == as.numeric(union3[i,1])))
  if(length(index) == 0){
    color1[i] = "blue"
  }
}

g3 = graph_from_data_frame(union3, directed = FALSE, vertices = subdata$ID)
E(g3)$color = color3
E(g3)$width = 2
V(g3)$color = lastname
V(g3)$shape[which(subdata$group=="A")] = "circle"
V(g3)$shape[which(subdata$group=="B")] = "square"
V(g3)$shape[which(subdata$group=="C")] = "rectangle"


#############################################################

par(mfrow=c(1,4))
plot(g, vertex.label=NA, layout = layout, vertex.size=7, main = "Predicted network \nwithout penalty adjustment")
plot(g1, vertex.label=NA, layout = layout, vertex.size=7, main = "Predicted network \nwithout penalty adjustment")
plot(g2, vertex.label=NA, layout = layout, vertex.size=7, main = "Predicted network \nwith alpha from greedy approach")
plot(g3, vertex.label=NA, layout = layout, vertex.size=7, main = "Predicted network \nwith alpha from Bayesian approach")
