library(readxl)

wages_2001 <- read_xls("data/bls_wages/oes01ma/MSA_2001_dl_1.xls") %>% select(area_name, occ_title, a_median) %>% mutate(year = 2001)
wages_2019 <- read_xlsx("data/bls_wages/oesm19ma/MSA_M2019_dl.xlsx") %>% select(area_name = area_title, occ_title, a_median) %>% mutate(year = 2019)

wages <-
  bind_rows(wages_2001, wages_2019) %>%
  filter(occ_title == "Industry Total")
  
area <- c("Sacramento, CA PMSA", "Los Angeles-Long Beach, CA PMSA")  

wages %>%
  filter(area_name == "Sacramento, CA PMSA")

wages_2001 %>% filter(area_name == "Los Angeles-Long Beach, CA PMSA") %>% filter(occ_title == "Industry Total") %>% select(a_wpct10, a_median, a_wpct90)
wages_2019 %>% filter(area_title == "Los Angeles-Long Beach-Anaheim, CA") %>% filter(occ_title == "All Occupations") %>% select(a_pct10, a_median, a_pct90)
