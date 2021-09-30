## load data
data("Bundesliga")
n <- nrow(Bundesliga)

## fit linear model
Bundesliga$logMarketValue <- log(Bundesliga$MarketValue)
fit <- lm(logMarketValue ~ Contract + Matches + Goals + Assists, 
          data=Bundesliga)

## perform K-fold cross-validation
perry(fit, foldControl(K = 5, R = 10), seed = 1234)

## perform random splitting
perry(fit, splitControl(m = n/3, R = 10), seed = 1234)

## perform bootstrap prediction error estimation
# 0.632 estimator
perry(fit, bootControl(R = 10, type = "0.632"), seed = 1234)
# out-of-bag estimator
perry(fit, bootControl(R = 10, type = "out-of-bag"), seed = 1234)
