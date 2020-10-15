library(tidyverse)
library(lubridate)
library(readxl)
library(stringr)
library(tidyquant)

prices3_br <-
  read_csv("data/zillow_historical_rents_and_prices/Metro_zhvi_bdrmcnt_3_uc_sfrcondo_tier_0.33_0.67_sm_sa_mon.csv") %>%
  pivot_longer(c(-RegionID, -SizeRank, -RegionName, -RegionType, -StateName),
               names_to = "date", values_to = "price_3br") %>%
  mutate(date = ymd(date)) %>%
  mutate(area = str_extract(RegionName, "[a-zA-Z \\-\\.]*"))

prices3_br %>%
  filter(SizeRank < 10) %>%
  ggplot(aes(date, price_3br / 1000, color = RegionName)) +
  geom_line()



rents_2001 <-
  read_xls("data/hud_historical_rents/FY2001_50th.xls") %>%
  mutate(type = str_extract(ID_AGIS2, "[a-zA-Z]*")) %>%
  filter(type == "MSA") %>%
  mutate(area = str_extract(AREANAME, "[a-zA-Z \\-\\.]*")) %>%
  select(state = ST, area,
         rent_0br = RENT0BR, rent_1br = RENT1BR, rent_2br = RENT2BR,
         rent_3br = RENT3BR, rent_4br = RENT4BR) %>%
  mutate(year = 2001)

rents_2002 <-
  read_xls("data/hud_historical_rents/FY2002_50th.xls") %>%
  mutate(area = str_extract(AREANAME, "[a-zA-Z \\-\\.]*")) %>%
  mutate(area = str_replace(area, "\\-\\-", "\\-")) %>%
  mutate(state = str_extract(AREANAME, "[A-Z][A-Z]")) %>%
  select(state = state, area,
         rent_0br = RENT_0, rent_1br = RENT_1, rent_2br = RENT_2,
         rent_3br = RENT_3, rent_4br = RENT_4) %>%
  mutate(year = 2002)

rents_2003 <-
  read_xls("data/hud_historical_rents/FY2003_50Rents_Area.xls") %>%
  mutate(area = str_extract(msaname, "[a-zA-Z \\-\\.]*")) %>%
  mutate(area = str_replace(area, "\\-\\-", "\\-")) %>%
  select(state = State_Alpha, area,
         rent_0br = rent50_0, rent_1br = rent50_1, rent_2br = rent50_2,
         rent_3br = rent50_3, rent_4br = rent50_4) %>%
  mutate(year = 2003)

rents_2004 <-
  read_xls("data/hud_historical_rents/FMR2004F50_FMRArea.xls") %>%
  mutate(area = str_extract(MSAName, "[a-zA-Z \\-\\.]*")) %>%
  mutate(area = str_replace(area, "\\-\\-", "\\-")) %>%
  select(state = State_Alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2004)

rents_2005 <-
  read_xls("data/hud_historical_rents/FMR_Area_50.xls") %>%
  mutate(area = str_extract(MSAName, "[a-zA-Z \\-\\.]*")) %>%
  mutate(area = str_replace(area, "\\-\\-", "\\-")) %>%
  select(state = State_Alpha, area,
         rent_0br = Rent50_0Bed, rent_1br = Rent50_1Bed, rent_2br = Rent50_2Bed,
         rent_3br = Rent50_3Bed, rent_4br = Rent50_4Bed) %>%
  mutate(year = 2005)

rents_2006 <-
  read_xls("data/hud_historical_rents/FY2006_Area_50th.xls") %>%
  mutate(area = str_extract(areaname, "[a-zA-Z \\-\\.]*")) %>%
  mutate(area = str_replace(area, "\\-\\-", "\\-")) %>%
  select(state = state_alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2006)

rents_2007 <-
  read_xls("data/hud_historical_rents/FY2007_Area_50th.xls") %>%
  mutate(area = str_extract(areaname, "[a-zA-Z \\-\\.]*")) %>%
  mutate(area = str_replace(area, "\\-\\-", "\\-")) %>%
  select(state = state_alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2007)

rents_2008 <-
  read_xls("data/hud_historical_rents/FY2008_area_50th_r.xls") %>%
  mutate(area = str_extract(Areaname, "[a-zA-Z \\-\\.]*")) %>%
  mutate(area = str_replace(area, "\\-\\-", "\\-")) %>%
  select(state = state_alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2008)

rents_2009 <-
  read_xls("data/hud_historical_rents/FMRArea_FY2009_50_Rev_Final.xls") %>%
  mutate(area = str_extract(Areaname, "[a-zA-Z \\-\\.]*")) %>%
  mutate(area = str_replace(area, "\\-\\-", "\\-")) %>%
  select(state = state_alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2009)

rents_2010 <-
  read_xls("data/hud_historical_rents/FMRArea_FY2010_50_Final.xls") %>%
  mutate(area = str_extract(Areaname, "[a-zA-Z \\-\\.]*")) %>%
  select(state = state_alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2010)

rents_2011 <-
  read_xls("data/hud_historical_rents/FMRArea_FY2011_50_Final.xls") %>%
  mutate(area = str_extract(Areaname, "[a-zA-Z \\-\\.]*")) %>%
  select(state = state_alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2011)

rents_2012 <-
  read_xls("data/hud_historical_rents/FY2012_FMRS_50_Area.xls") %>%
  mutate(area = str_extract(Areaname, "[a-zA-Z \\-\\.]*")) %>%
  select(state = state_alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2012)

rents_2013 <-
  read_xls("data/hud_historical_rents/FMRArea_FY2013_50_Final.xls") %>%
  mutate(area = str_extract(Areaname, "[a-zA-Z \\-\\.]*")) %>%
  select(state = state_alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2013)

rents_2014 <-
  read_xls("data/hud_historical_rents/FMRArea_FY2014_50_RevFinal.xls") %>%
  mutate(area = str_extract(Areaname, "[a-zA-Z \\-\\.]*")) %>%
  select(state = state_alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2014)

rents_2015 <-
  read_xls("data/hud_historical_rents/FMRArea_FY2015_50_RevFinal.xls") %>%
  mutate(area = str_extract(Areaname, "[a-zA-Z \\-\\.]*")) %>%
  select(state = state_alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2015)

rents_2016 <-
  read_xlsx("data/hud_historical_rents/FMRArea_FY2016F_50_RevFinal.xlsx") %>%
  mutate(area = str_extract(areaname, "[a-zA-Z \\-\\.]*")) %>%
  select(state = state_alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2016)

rents_2017 <-
  read_xlsx("data/hud_historical_rents/FMRArea_FY2017_50_rev.xlsx") %>%
  mutate(area = str_extract(areaname, "[a-zA-Z \\-\\.]*")) %>%
  select(state = state_alpha, area,
         rent_0br = Rent50_0, rent_1br = Rent50_1, rent_2br = Rent50_2,
         rent_3br = Rent50_3, rent_4br = Rent50_4) %>%
  mutate(year = 2017)

rents_2018 <-
  read_xlsx("data/hud_historical_rents/FY2018_50_FMRArea_rev.xlsx") %>%
  mutate(area = str_extract(areaname18, "[a-zA-Z \\-\\.]*")) %>%
  select(state = state_alpha, area,
         rent_0br = rent50_0, rent_1br = rent50_1, rent_2br = rent50_2,
         rent_3br = rent50_3, rent_4br = rent50_4) %>%
  mutate(year = 2018)

rents_2019 <-
  read_xlsx("data/hud_historical_rents/FY2019_50_FMRArea_rev.xlsx") %>%
  mutate(area = str_extract(areaname19, "[a-zA-Z \\-\\.]*")) %>%
  select(state = state_alpha, area,
         rent_0br = rent50_0, rent_1br = rent50_1, rent_2br = rent50_2,
         rent_3br = rent50_3, rent_4br = rent50_4) %>%
  mutate(year = 2019)

rents_2020 <-
  read_xlsx("data/hud_historical_rents/FY2020_50_FMRArea_rev2.xlsx") %>%
  mutate(area = str_extract(areaname20, "[a-zA-Z \\-\\.]*")) %>%
  select(state = state_alpha, area,
         rent_0br = rent50_0, rent_1br = rent50_1, rent_2br = rent50_2,
         rent_3br = rent50_3, rent_4br = rent50_4) %>%
  mutate(year = 2020)

rents <- bind_rows(rents_2001, rents_2002, rents_2003, rents_2004,
                   rents_2005, rents_2006, rents_2007, rents_2008, rents_2009,
                   rents_2010, rents_2011, rents_2012, rents_2013, rents_2014,
                   rents_2015, rents_2016, rents_2017, rents_2018, rents_2019,
                   rents_2020)

rents %>%
  group_by(area) %>%
  filter(n() == 5) %>%
  ungroup() %>%
  # filter(area %in% sample(unique(.$area), 9)) %>%
  filter(state %in% c("CA", "TN", "AZ", "OR", "WA", "NY", "FL", "IN", "TX")) %>%
  # mutate(area = paste(area, state)) %>%
  ggplot(aes(year, rent_3br, group = area)) +
  geom_point() +
  geom_line() +
  facet_wrap(state ~ .)

rents %>%
  filter(year %in% c(2010, 2020)) %>%
  arrange(state, area, year) %>%
  # distinct(state, area, .keep_all = TRUE) %>%
  pivot_wider(id_cols = c(state, area), names_from = year, values_from = rent_3br) %>%
  filter(!is.na(`2010`) & !is.na(`2020`)) %>%
  mutate(change = `2020` / `2010`) %>%
  # filter(state %in% c("CA", "TN", "AZ", "OR", "WA", "NY", "FL", "IN", "TX")) %>%
  ggplot(aes(change)) +
  geom_histogram() +
  facet_wrap(state ~ .)

rents %>%
  filter(state == "CA") %>%
  select(state, area, year, rent_3br) %>%
  ggplot(aes(year, rent_3br, color = area)) +
  geom_point() +
  geom_line() +
  facet_wrap(area ~ .)

# Compare rents and prices
prices_2020 <- prices3_br %>% filter(date == "2020-06-30") %>% select(state = StateName, area, price_3br)
rents_2020 <- rents %>% filter(year == 2020) %>% select(state, area, rent_3br)

i <- 0.03 / 12
n <- 360
inner_join(prices_2020, rents_2020) %>%
  mutate(principal = price_3br * 0.9, downpayment = price_3br - principal) %>%
  mutate(mortgage = principal * (i * (1+i)^n) / ((1+i)^n - 1)) %>%
  mutate(profit = rent_3br - mortgage) %>%
  arrange(desc(profit)) %>%
  ggplot(aes(profit)) +
  geom_histogram()

# Redfin
redfin <- read_tsv("data/redfin/weekly_housing_market_data_most_recent.tsv")

# Zillow sale prices
sale_prices <- read_csv("data/zillow_historical_rents_and_prices/Metro_median_sale_price_uc_SFRCondo_raw_month.csv")

sale_prices %>%
  filter(SizeRank < 5) %>%
  select(-RegionID, -SizeRank, -RegionType, -StateName) %>%
  pivot_longer(-RegionName, names_to = "date", values_to = "med_sale_price")

# 3br vs stocks
sp500 <-
  tq_get("SPY", get = "stock.prices", from = "2001-01-01", to = "2020-02-01") %>%
  filter(day(date) == 2) %>%
  select(date, price = close) %>%
  mutate(asset = "SPY") %>%
  mutate(price = price) %>%
  mutate(rent_or_buy = "buy")

price_3br <-
  prices3_br %>%
  filter(RegionName %in% c("Sacramento, CA", "Los Angeles-Long Beach-Anaheim, CA", "Austin, TX", "Toledo, OH")) %>%
  filter(date >= "2001-01-01") %>%
  inner_join(tibble(RegionName = c("Sacramento, CA", "Los Angeles-Long Beach-Anaheim, CA", "Austin, TX", "Toledo, OH"),
                    common_name = c("Sacramento", "LA-LB", "Austin", "Toledo"))) %>%
  mutate(asset = paste("Buy", common_name)) %>%
  select(date, price = price_3br, asset = common_name) %>%
  mutate(rent_or_buy = "buy")

rent_3br <-
  rents %>%
  filter(area %in% c("Sacramento", "Sacramento--Arden-Arcade--Roseville", "Sacramento--Roseville--Arden-Arcade",
                     "Los Angeles-Long Beach", "Los Angeles-Long Beach-Glendale",
                     "Austin-San Marcos", "Austin-Round Rock-San Marcos", "Austin-Round Rock",
                     "Toledo")) %>%
  inner_join(tibble(area = c("Sacramento", "Sacramento--Arden-Arcade--Roseville", "Sacramento--Roseville--Arden-Arcade",
                             "Los Angeles-Long Beach", "Los Angeles-Long Beach-Glendale",
                             "Austin-San Marcos", "Austin-Round Rock-San Marcos", "Austin-Round Rock",
                             "Toledo"),
                    common_name = c(rep("Sacramento", 3), rep("LA-LB", 2), rep("Austin", 3), rep("Toledo",1)))) %>%
  mutate(date = ymd(paste0(year, "-01-01"))) %>%
  select(date, price = rent_3br, asset = common_name) %>%
  mutate(rent_or_buy = "rent")

bind_rows(sp500, price_3br, rent_3br) %>%
  filter(date >= "2015-01-01") %>%
  arrange(date) %>%
  group_by(asset, rent_or_buy) %>%
  mutate(price = price / head(price, 1)) %>%
  ungroup() %>%
  ggplot(aes(date, price, color = asset)) +
  geom_line(aes(linetype = rent_or_buy))


#### Cap rate for various cities
rent_3br <-
  rents %>%
  filter(area %in% c("Sacramento", "Sacramento--Arden-Arcade--Roseville", "Sacramento--Roseville--Arden-Arcade",
                     "Los Angeles-Long Beach", "Los Angeles-Long Beach-Glendale",
                     "Austin-San Marcos", "Austin-Round Rock-San Marcos", "Austin-Round Rock",
                     "Toledo")) %>%
  inner_join(tibble(area = c("Sacramento", "Sacramento--Arden-Arcade--Roseville", "Sacramento--Roseville--Arden-Arcade",
                             "Los Angeles-Long Beach", "Los Angeles-Long Beach-Glendale",
                             "Austin-San Marcos", "Austin-Round Rock-San Marcos", "Austin-Round Rock",
                             "Toledo"),
                    common_name = c(rep("Sacramento", 3), rep("LA-LB", 2), rep("Austin", 3), rep("Toledo",1)))) %>%
  mutate(date = ymd(paste0(year, "-01-01"))) %>%
  select(year, monthly_rent = rent_3br, common_name)

price_3br <-
  prices3_br %>%
  filter(RegionName %in% c("Sacramento, CA", "Los Angeles-Long Beach-Anaheim, CA", "Austin, TX", "Toledo, OH")) %>%
  filter(date >= "2001-01-01") %>%
  inner_join(tibble(RegionName = c("Sacramento, CA", "Los Angeles-Long Beach-Anaheim, CA", "Austin, TX", "Toledo, OH"),
                    common_name = c("Sacramento", "LA-LB", "Austin", "Toledo"))) %>%
  mutate(asset = paste("Buy", common_name)) %>%
  filter(month(date) == 1, day(date) == 31) %>%
  mutate(year = year(date)) %>%
  select(year, price = price_3br, common_name)

mortgage_rates <-
  read_csv("data/MORTGAGE30US.csv") %>%
  filter(month(DATE) == 1) %>%
  mutate(year = year(DATE)) %>%
  group_by(year) %>%
  filter(DATE == min(DATE)) %>%
  select(year, mortgage_rate_30_year = MORTGAGE30US) %>%
  ungroup()

median_us_home <-
  read_csv("data/MSPUS.csv") %>%
  filter(month(DATE) == 1) %>%
  mutate(year = year(DATE)) %>%
  group_by(year) %>%
  filter(DATE == min(DATE)) %>%
  select(year, med_house_price = MSPUS) %>%
  ungroup()

median_us_home %>%
  inner_join(mortgage_rates) %>%
  mutate(i = mortgage_rate_30_year / 1200) %>%
  mutate(monthly_mortgage = (med_house_price - 20000) * (i * (1+i)^n) / ((1+i)^n - 1)) %>%
  select(year, med_house_price, mortgage_rate_30_year, monthly_mortgage) %>%
  pivot_longer(-year) %>%
  ggplot(aes(year, value)) +
  geom_line() +
  facet_grid(name ~ ., scales = "free")
  

n <- 360
rent_mortgage <-
  inner_join(rent_3br, price_3br) %>%
  inner_join(mortgage_rates) %>%
  mutate(principal = price * 0.8,
         downpayment = price - principal,
         i = (mortgage_rate_30_year / 100) / 12) %>%
  mutate(mortgage = principal * (i * (1+i)^n) / ((1+i)^n - 1)) %>%
  mutate(profit =  monthly_rent - mortgage)

rent_mortgage %>%
  ggplot(aes(year, profit, color = common_name)) +
  geom_point() +
  geom_line()
  

# Rent vs buying
1.39e6*0.9 * (i * (1+i)^n) / ((1+i)^n - 1)

df <-
  tibble(month = c(1:360, 361:(361+179))) %>%
  mutate(mortgage = ifelse(month <= 360, 5274, 0),
         property_tax = 750,
         insurance = 100,
         maintenance = 1000) %>%
  mutate(buy_total = mortgage + property_tax + insurance + maintenance) %>%
  mutate(rent = 4545 * 1.003357^(month - 1)) %>%
  #mutate(rent = 4545) %>%
  mutate(rent_savings = buy_total - rent) %>%
  mutate(house_savings = rent - buy_total) %>%
  mutate(rent_savings = ifelse(rent_savings < 0, 0, rent_savings)) %>%
  mutate(house_savings = ifelse(house_savings < 0, 0, house_savings)) %>%
  mutate(rent_savings_pv = rent_savings * 1.007974^(540 - month)) %>%
  mutate(house_savings_pv = house_savings * 1.007974^(540 - month))

sum(df$rent_savings_pv) + (1.39e6*0.15)*1.1^45 # RENTING
1.39e6 * 1.064139^45 + sum(df$house_savings_pv) # BUYING

# Angad culver city home
cash_flow_analysis <- function(house_price, realtor_closing_pct, down_payment_pct,
                               monthly_house_apprec_rate, monthly_stock_apprec_rate, yearly_property_tax,
                               yearly_insurance, yearly_maintenance, n = 360, i = 0.03 / 12,
                               rent, monthly_rent_appreciation, vacancy_rate) {
  
  down_payment <- down_payment_pct * house_price
  realtor_closing_cost <- house_price * realtor_closing_pct
  total_upfront_cost <- down_payment + realtor_closing_cost
  
  loan_amount <- (1 - down_payment_pct) * house_price
  monthly_mortgage <- loan_amount * (i * (1+i)^n) / ((1+i)^n - 1)
  
  df <-
    tibble(month = 1:n) %>%
    mutate(mortgage = monthly_mortgage) %>%
    mutate(principal_part = (mortgage - i * loan_amount) * (1 + i)^(month - 1)) %>%
    mutate(total_princ_paid = cumsum(principal_part)) %>%
    mutate(loan_balance = pmax(loan_amount - total_princ_paid, 0)) %>%
    mutate(home_value = house_price * (1 + monthly_house_apprec_rate)^month) %>%
    mutate(total_equity = home_value - loan_balance) %>%
    mutate(stock_value = total_upfront_cost * (1 + monthly_stock_apprec_rate)^month) %>%
    mutate(property_tax = yearly_property_tax / 12,
           insurance = yearly_insurance / 12,
           maintenance = yearly_maintenance / 12) %>%
    mutate(monthly_cost = mortgage + property_tax + insurance + maintenance) %>%
    mutate(rent = rent * (1 + monthly_rent_appreciation)^(month - 1) * (1 - vacancy_rate)) %>%
    mutate(cash_flow = rent - monthly_cost) %>%
    mutate(total_earnings = cumsum(cash_flow)) %>%
    mutate(earnings_plus_equity = total_equity + total_earnings) %>%
    mutate()
  
  list(df = df,
       down_payment = down_payment,
       realtor_closing_cost = realtor_closing_cost,
       total_upfront_cost = total_upfront_cost,
       monthly_mortgage = monthly_mortgage,
       monthly_cost = df$monthly_cost[1],
       monthly_cash_flow = df$cash_flow[1])
}

angad_culver_city <-
  cash_flow_analysis(house_price = 1.9e6,
                   realtor_closing_pct = 0.05,
                   down_payment_pct = 0.2,
                   monthly_house_apprec_rate = 0.0046,
                   monthly_stock_apprec_rate = 0.008,
                   yearly_property_tax = 14345,
                   yearly_insurance = 2000,
                   yearly_maintenance = 10000,
                   n = 360,
                   i = 0.03 / 12,
                   rent = 8342,
                   monthly_rent_appreciation = 0.0016,
                   vacancy_rate = 0.05)

angad_culver_city$df %>%
  select(month, stock_value, total_equity, earnings_plus_equity, total_earnings) %>%
  mutate(years = month / 12) %>%
  select(-month) %>%
  pivot_longer(-years) %>%
  ggplot(aes(years, value, color = name)) +
  geom_line()

# West street duplex
west_street <-
  cash_flow_analysis(house_price = 4.89e5,
                     realtor_closing_pct = 0.05,
                     down_payment_pct = 0.2,
                     monthly_house_apprec_rate = 0.0046,
                     monthly_stock_apprec_rate = 0.008,
                     yearly_property_tax = 14345,
                     yearly_insurance = 2000,
                     yearly_maintenance = 10000,
                     n = 360,
                     i = 0.03 / 12,
                     rent = 8342,
                     monthly_rent_appreciation = 0.0016,
                     vacancy_rate = 0.05)