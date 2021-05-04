
library(dplyr)

path <- "/SFS/scratch/zhanyilo/actg_rtb_50/"
result = result_true = list()

for(i in 1:1000){
  try({
    load(file.path(path, paste0(i, ".Rdata")))

    result[[i]] <- res
    result_true[[i]] <- res_true

  })
}

result <- bind_rows(result)
result_true <- bind_rows(result_true)

s_true <- result_true %>% group_by(group) %>% summarise_all(mean) %>% select(-n, - n_mi)

res1 <- result %>% left_join(s_true) %>% subset(delta != 88) %>%

  group_by(n_mi, group, delta) %>%
  mutate(true = mean(true)) %>% ungroup()  %>%
  mutate(
    covr    = (true > rmst - 1.96 * sd) & (true < rmst + 1.96 * sd),
    covr_wb = (true > rmst - 1.96 * wb_sd) & (true < rmst + 1.96 * wb_sd)) %>%
  group_by(n, n_mi, group, delta) %>%
  summarise(est = mean(rmst),
            true = mean(true),
            sd_true = sd(rmst),
            sd_est = mean(sd),
            wb_sd_est = mean(wb_sd),
            covr = mean(covr),
            covr_wb = mean(covr_wb))

res1 %>% subset(n_mi == 50 & n == 1000 & delta == 2)
res1 %>% subset(n_mi == 20 & n == 1000 & delta == 2)

res2 <- result %>% left_join(s_true) %>% subset(delta == 88) %>%

  group_by(n, group, delta) %>%
  mutate(true = mean(true)) %>% ungroup()  %>%
  mutate(
    covr    = (true > rmst - 1.96 * sd) & (true < rmst + 1.96 * sd),
    covr_wb = (true > rmst - 1.96 * wb_sd) & (true < rmst + 1.96 * wb_sd)) %>%
  group_by(n, group, delta) %>%
  summarise(est = mean(rmst),
            true = mean(true),
            sd_true = sd(rmst),
            wb_sd_est = mean(sd),
            covr_wb =  mean(covr))


res <- res1 %>% subset(n_mi %in% c(10, 20, 50) | delta == 88)
res3 <- res
count(ungroup(res3), delta)
count(ungroup(res), delta)

# res <- bind_rows(res1, res2) %>% subset(n_mi %in% c(20) | delta == 88)
res2
res4 <- res3 %>% subset(! (n == 1000 & n_mi == 50))
res5 <- subset(res1, n== 1000 & n_mi == 50)
res <- bind_rows(res4, res2, res5)
write.csv(res, file = "simulation/atcg.csv")
