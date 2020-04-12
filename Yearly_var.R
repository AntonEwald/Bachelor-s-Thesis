#------------------- Code to get the lowest and highest yearly variances -----------------
## @Author : Anton Holm

target_locations <- metals %>% group_by(LOC, YEAR) %>% count() %>% filter(n>9) %>% ungroup() %>% select(LOC) %>% filter(LOC != "Utlängan (spring)") %>% filter(LOC !="Ängskärsklubb (spring)") %>% filter(LOC != "Örefjärden")


locations_w_many_obs <- metals %>% filter(LOC %in% target_locations[[1]])

min_std_yearly <- locations_w_many_obs %>% group_by(LOC, YEAR) %>%  summarise(std = sd(NI)) %>% na.omit() %>% ungroup() %>% group_by(YEAR) %>% summarise(min_sd = min(std))

max_std_yearly <- locations_w_many_obs %>% group_by(LOC, YEAR) %>%  summarise(std = sd(NI)) %>% na.omit() %>% ungroup() %>% group_by(YEAR) %>% summarise(max_sd = max(std))
yearly_variances <- cbind(min_std_yearly[,2], max_std_yearly[,2])


save(yearly_variances, file="yearly_variances.Rda")