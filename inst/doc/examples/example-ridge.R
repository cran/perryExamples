## load data
data("Bundesliga")
Bundesliga <- Bundesliga[, -(1:2)]
f <- log(MarketValue) ~ Age + I(Age^2) + .
mf <- model.frame(f, data=Bundesliga)
x <- model.matrix(terms(mf), mf)[, -1]
y <- model.response(mf)

## set up repeated random splits
splits <- splitControl(m = 40, R = 10)

## select optimal penalty parameter
lambda <- seq(600, 0, length.out = 50)
fit <- ridge(x, y, lambda = lambda, splits = splits, seed = 2014)
fit

## plot prediction error results
plot(fit, method = "line")
