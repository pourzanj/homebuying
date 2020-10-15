library(tidyverse)

get_zipcode_results <- function(zipcode) {
  Sys.sleep(exp(rnorm(1)))
  craigr::rentals(location = "sfbay", area = "eby",
                  max_results = 400,
                  postal = zipcode, search_distance = 2.0,
                  bedrooms = 1, posted_today = FALSE) %>%
    mutate(postal = zipcode)
}

zipcodes <-
  tribble(
    ~station,               ~zipcode,
    "West Oakland",         94607,
    "MacArthur",            94609,
    "Rockridge",            94618,
    "Downtown Berkeley",    94704,
    "Ashby",                94703,
    "Orinda",               94563,
    "El Cerrito",           94530,
    "El Cerrito Del Norte", 94530,
    "North Berkeley",       94702,
    "Lafayette",            94549
)

results <-
  zipcodes$zipcode %>%
  map(get_zipcode_results) %>%
  bind_rows()

# write_csv(results, "data/historical_prices.csv")

results %>%
  inner_join(zipcodes, by = c("postal" = "zipcode")) %>%
  ggplot(aes(Price)) +
  geom_histogram() +
  facet_wrap(~station)
