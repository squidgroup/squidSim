### testing known_predictors

#devtools::install("./github/squidSim")
library(squidSim)

k_p <- MASS::mvrnorm(100,c(0,0),diag(2))
length(colnames(k_p))
colnames(k_p) <- c("x1","x2")

squid_dat <- simulate_population(
	n=100,
	parameters = list(
    observation =list(
      names = c("temperature","rainfall"),
      beta = c(0.5,0.3)
    ),
    residual = list(
      vcov = 0.3
    )
  ),
	known_predictors=list(predictors=cbind(x1=k_p[,1]), beta=c(x2=2))
	)


head(squid_dat$pred[[1]])

squid_dat$known
dat <- get_population_data(squid_dat)

lm(y~x1+x2,dat)




data_structure=make_structure("sex(2)",repeat_obs=100, 
### testing factors

squid_dat <- simulate_population(
	data_structure=make_structure("sex(2)",repeat_obs=100, level_names=list(sex=c("Male","Female"))),
	parameters = list(
    sex =list(
    	fixed=TRUE,
      names = c("Male","Female"),
      beta = c(0.5,0.3)
    ),
    residual = list(
      vcov = 0.3
    )
  )
)
dat <- get_population_data(squid_dat)
head(dat)
boxplot(y~sex,dat)

