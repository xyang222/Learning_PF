load("data_community_large.Rdata")

combination = combn(nrow(data), 2)
npeople = nrow(data)
#################################################birth/death
compare_date = function(people_data, combination){
  birth1 = people_data$birth[combination[1,]]
  death1 = people_data$death[combination[1,]]
  birth2 = people_data$birth[combination[2,]]
  death2 = people_data$death[combination[2,]]
  return(birth1 < death2 & death1 > birth2)
}

overlap_l = compare_date(data, combination)
overlap_m = matrix(rep(TRUE, npeople*npeople), nrow = npeople)
overlap_m[lower.tri(overlap_m)] = overlap_l
overlap_m = t(overlap_m)
overlap_m[lower.tri(overlap_m)] = overlap_l

#################################################lastname

compare_lastname = function(people_data, combination){
  lastname1 = people_data$last_names[combination[1,]]
  lastname2 = people_data$last_names[combination[2,]]
  return(lastname1 == lastname2)
}

family_l = compare_lastname(data, combination)
family_m = matrix(rep(TRUE, npeople*npeople), nrow = npeople)
family_m[lower.tri(family_m)] = family_l
family_m = t(family_m)
family_m[lower.tri(family_m)] = family_l

#################################################group
compare_group = function(ingroup, combination){
  group1 = ingroup[combination[1,]]
  group2 = ingroup[combination[2,]]
  return(group1 & group2)
}
group = matrix(rep(FALSE, npeople*npeople), nrow = npeople)

add_group = function(group_id, group){
  ingroup = (data$group == group_id)

  group_l = compare_group(ingroup, combination)
  group_m = matrix(rep(TRUE, npeople*npeople), nrow = npeople)
  group_m[lower.tri(group_m)] = group_l
  group_m = t(group_m)
  group_m[lower.tri(group_m)] = group_l

  group = (group | group_m)
  return(group)
}

group_A = add_group("A", group)
group_B = add_group("B", group)
group_C = add_group("C", group)
###########################################################

save(overlap_m, family_m, group_A, group_B, group_C, file = "penalties_large.Rdata")
