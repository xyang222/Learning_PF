library(gender)

load("sampmatrix1.Rdata")
load("ODNB_dataset.Rdata")

IDs = which(nodeset$ext_birth >= 1500 & nodeset$ext_death <= 1575)
people = nodeset[IDs,]
doc_counts = dcmat[,IDs]
doc_counts = doc_counts[which(rowSums(doc_counts) != 0),]
emptypeople = which(colSums(doc_counts) != 0)
people = people[emptypeople, ]
doc_counts = doc_counts[,emptypeople]
dim(doc_counts)

relationship = read.csv("SDFB_relationships_2019_03_08.csv")
relationship = relationship[which(relationship$person1_index %in% people$SDFB_ID & relationship$person2_index %in% people$SDFB_ID),]
relationship = relationship[which(relationship$original_certainty >= 50),]
relationship = relationship[which(relationship$created_by == 2),]
##########################################################
#table(people$occupation)

occupation = tolower(people$occupation)
occupation = strsplit(occupation, " ")
poetOrNot = lapply(occupation, function(list){
  if("poet" %in% list | "writer" %in% list | "author" %in% list)
    return(TRUE)
  else
    return(FALSE)
  })
poetOrNot = unlist(poetOrNot)

churchOrNot = lapply(occupation, function(list){
  if("church" %in% list | "religious" %in% list |
     "bishop" %in% list | "catholic" %in% list)
    return(TRUE)
  else
    return(FALSE)
})
churchOrNot = unlist(churchOrNot)

royalOrNot = lapply(occupation, function(list){
  if("royal" %in% list | "king" %in% list |
     "queen" %in% list | "regent" %in% list)
    return(TRUE)
  else
    return(FALSE)
})
royalOrNot = unlist(royalOrNot)

# predictedGender = sapply(people$first_name, function(name){
#  return(gender(name, countries = "United Kingdom", method = "napp")$gender)})
# predictedGender = unlist(predictedGender)
# maleOrNot = (predictedGender == "male")

###############################################################
combination = combn(nrow(people), 2)
npeople = nrow(people)

compare_date = function(people_data, combination){
  birth1 = people_data$ext_birth[combination[1,]]
  death1 = people_data$ext_death[combination[1,]]
  birth2 = people_data$ext_birth[combination[2,]]
  death2 = people_data$ext_death[combination[2,]]
  return(birth1 < death2 & death1 > birth2)
}
overlap_l = compare_date(people, combination)
overlap_m = matrix(rep(TRUE, npeople*npeople), nrow = npeople)
overlap_m[lower.tri(overlap_m)] = overlap_l
overlap_m = t(overlap_m)
overlap_m[lower.tri(overlap_m)] = overlap_l

compare_lastname = function(people_data, combination){
  lastname1 = people_data$surname[combination[1,]]
  lastname2 = people_data$surname[combination[2,]]
  return(lastname1 == lastname2)
}
family_l = compare_lastname(people, combination)
family_m = matrix(rep(TRUE, npeople*npeople), nrow = npeople)
family_m[lower.tri(family_m)] = family_l
family_m = t(family_m)
family_m[lower.tri(family_m)] = family_l

compare_group = function(ingroup, combination){
  group1 = ingroup[combination[1,]]
  group2 = ingroup[combination[2,]]
  return(group1 & group2)
}
group = matrix(rep(FALSE, npeople*npeople), nrow = npeople)

add_group = function(covariate, group){
  ingroup = covariate

  group_l = compare_group(ingroup, combination)
  group_m = matrix(rep(TRUE, npeople*npeople), nrow = npeople)
  group_m[lower.tri(group_m)] = group_l
  group_m = t(group_m)
  group_m[lower.tri(group_m)] = group_l

  group = (group | group_m)
  return(group)
}

group_poet = add_group(poetOrNot, group)
group_church = add_group(churchOrNot, group)
group_royal = add_group(royalOrNot, group)

#group_male = add_group(maleOrNot, group)
#group_female = 1-group_male

save(people, relationship, doc_counts, overlap_m, family_m, group_poet, group_church, group_royal, file = "realData.Rdata")
