library(randomNames)
library(igraph)

compare_date = function(people_data, combination){
  birth1 = people_data$birth[combination[1,]]
  death1 = people_data$death[combination[1,]]
  birth2 = people_data$birth[combination[2,]]
  death2 = people_data$death[combination[2,]]
  return(birth1 < death2 & death1 > birth2)
}

###################################################

set.seed(1)
lastName = randomNames(30, ethnicity = rep(5, 20), which.names = "last")

data = data.frame()
relationship = data.frame()

set.seed(83)
for(i in 1:50){
  #last = lastName[i]
  last = sample(lastName, 1)
  size = sample(3:15, 1)
  firstnames = randomNames(size, ethnicity = rep(5, size), which.names = "first", return.complete.data = TRUE)
  birth = sample(1500:1600, size)
  death = birth+sample(5:70, size)
  family = data.frame(ID = i*100 + (1:size) , first_names = firstnames[,3], last_names = rep(last, size), gender = firstnames[,1], birth = birth, death = death)
  data = rbind(data, family)

  combination = combn(1:size,2)
  overlap = compare_date(family, combination)
  relation = matrix(family$ID[combination[,which(overlap)]], nrow = 2)
  relation_sample = sample(1:ncol(relation), ceiling(ncol(relation)*0.5))
  relationship = rbind(relationship, t(relation[,relation_sample]))
}

data$group = sample(c("A", "B", "C"), nrow(data), replace = TRUE, prob = c(0.5, 0.25, 0.25))

combination = combn(1:nrow(data),2)
overlap = compare_date(data, combination)
potential_relation = matrix(data$ID[combination[,which(overlap)]], nrow = 2)

relation_all = potential_relation[,sample(1:ncol(potential_relation), 300)]

relation_A_group = potential_relation[,which(data$group[match(potential_relation[1,], data$ID)]== "A" & data$group[match(potential_relation[2,], data$ID)]== "A")]
relation_A = relation_A_group[,sample(1:ncol(relation_A_group), 100)]

relation_B_group = potential_relation[,which(data$group[match(potential_relation[1,], data$ID)]== "B" & data$group[match(potential_relation[2,], data$ID)]== "B")]
relation_B = relation_B_group[,sample(1:ncol(relation_B_group), 100)]

relation_C_group = potential_relation[,which(data$group[match(potential_relation[1,], data$ID)]== "C" & data$group[match(potential_relation[2,], data$ID)]== "C")]
relation_C = relation_C_group[,sample(1:ncol(relation_C_group), 50)]

relationship = rbind(relationship, t(relation_all), t(relation_A), t(relation_B),t(relation_C))
relationship = unique(relationship)
relationship = relationship[!duplicated(apply(relationship,1,function(x) paste(sort(x),collapse=''))),]

save(data, relationship, file = "data_community_large.Rdata")
