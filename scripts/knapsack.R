library(nilde)

# Function we're going to copy. a is cost of each item. n is max budget.
# it doesn't allow for non-integer costs
get.knapsack(c(0,0,0), a = c(1, 3, 4), n = 5)

get_all_possible_moves <- function (a, n) 
{
  if (length(a) < 2) {
    stop("length of vector 'a' has to be more than 1")
  }

  ra <- rank(a, ties.method = "first")
  a <- sort(a)
  l <- length(a)
  out <- numeric(0)
  a1 <- c(a[ra], 1)
  ra <- rank(a1, ties.method = "first")
  M <- floor(n/min(a1))
  a1 <- sort(a1)
  l1 <- length(a1)
  b <- c(floor(n/a1[l1]), rep(NA, l1 - 2))
  
  # Recursively list all solutions
  out <- recursive.list_knapsack(numeric(0), b, a1, n, M)
  
  # Organize back into matrix
  if (length(out) == 0) {
    out <- NULL
  }
  else {
    dim(out) <- c(l1, length(out)/l1)
    #out <- as.matrix(out[ra, ], l1, length(out)/l1)
  }
  
  out
}

# w keeps track of how many of each item you've already bought
# b keeps track of how many items are left and the highest quantity of the
# most expensive quantity you can buy
# a is how much each item costs
# n is the budget
recursive.list_knapsack <- function (w, b, a, n, M) 
{
  S <- rep(0, 0)
  d <- 0
  
  if (length(b)) {
    # You can buy up to b[1] of the current item. Iterate through all the
    # possibilities with 0, 1, 2, ..., b[1] of the current item and append
    for (i in seq(0, b[1])) {
      # If there are still items left figure out how much of the next highest
      # priced item you can buy if you buy i of the current item
      if (length(b) > 1) {
        d <- sum(w * a[length(a):(length(a) - length(w) + 
                                    1)]) + i * a[length(b) + 1]
        b[2] <- floor((n - d)/a[length(b)])
      }
      S <- c(S, Recall(c(w, i), b[-1], a, n, M))
    }
  }
  else {
    return(list_knapsack(rev(w), a, n, M))
  }
  S
}

list_knapsack <- function (ss, a, n, M) 
{
  res <- numeric(0)
  
  # left over cash
  s1 <- (n - sum(ss * a[-1]))
  #if (s1 >= 0 && s1 == floor(s1) && sum(c(s1, ss)) <= M) 
  res <- c(s1, ss)
  res
}

debugonce(get_all_possible_moves)
get_all_possible_moves(a = c(1.5, 3.2, 4.1), n = 5)

get_all_possible_moves(a = c(0.65, 0.7, 1.4, 2.35, , 4.1, 4.7), n = 5) %>% ncol()
