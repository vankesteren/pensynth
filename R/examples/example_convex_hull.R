# create some data
set.seed(45)
X0 <- matrix(runif(20), nrow = 2)
X1 <- matrix(c(.5, .5))

# test if X1 is in the convex hull:
in_convex_hull(X1, X0)

# also works with multiple units in X1:
X1 <- cbind(X1, c(1.3, -3))
in_convex_hull(X1, X0)
