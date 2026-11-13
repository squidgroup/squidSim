res <- list(residual = list(vcov = 1))
ped <- data.frame(id = as.character(1:6), dam = c(NA, NA, "1", "1", "3", "3"), sire = c(NA, NA, "2", "2", "4", "4"))

sim_col <- function(s, col) do.call(cbind, lapply(s$predictors, function(x) x[, col]))

test_that("repeated calls without a seed give different data", {
  set.seed(1)
  v <- replicate(500, simulate_population(n = 5, parameters = res)$y[[1]][1])
  expect_equal(length(unique(v)), 500)
})

test_that("a seed reproduces data and leaves the RNG stream alone", {
  a <- simulate_population(n = 5, parameters = res, seed = 42)
  b <- simulate_population(n = 5, parameters = res, seed = 42)
  expect_identical(a$y, b$y)
  set.seed(1); x <- runif(1)
  set.seed(1); simulate_population(n = 5, parameters = res, seed = 42); y <- runif(1)
  expect_identical(x, y)
})

test_that("sqrt link squares the linear predictor", {
  s <- simulate_population(n = 5000, family = "poisson", link = "sqrt", seed = 1,
    parameters = list(intercept = 3, residual = list(vcov = 1e-4)))
  expect_equal(mean(s$y[[1]]), 9, tolerance = 0.05)
})

test_that("zero variance predictors equal their mean", {
  s <- simulate_population(n = 100, seed = 1, parameters = list(
    observation = list(names = c("a", "b"), mean = c(5, 1), vcov = c(0, 1), beta = c(1, 1)),
    residual = list(vcov = 1)))
  expect_true(all(get_population_data(s)$a == 5))
})

test_that("singular vcov and 1x1 vcorr work", {
  s <- simulate_population(n = 100, seed = 1, parameters = list(
    observation = list(names = c("a", "b"), vcov = matrix(1, 2, 2), beta = c(1, 1)),
    residual = list(vcov = 1)))
  d <- get_population_data(s)
  expect_equal(cor(d$a, d$b), 1)
  s <- simulate_population(n = 5000, seed = 1, parameters = list(
    observation = list(names = "a", vcorr = matrix(2.25), beta = 1), residual = list(vcov = 1)))
  expect_equal(var(get_population_data(s)$a), 2.25, tolerance = 0.1)
})


test_that("mean is respected with pedigree, phylogeny and cov_str", {
  p <- function(nm) setNames(list(list(names = "g", mean = 5, vcov = 1, beta = 1), list(vcov = 1)), c(nm, "residual"))
  s <- simulate_population(data_structure = data.frame(animal = ped$id), pedigree = list(animal = ped),
    n_pop = 300, seed = 1, parameters = p("animal"))
  expect_equal(mean(sim_col(s, "g")), 5, tolerance = 0.05)

  # skip_if_not_installed("ape")
  # tr <- ape::rcoal(6)
  # s <- simulate_population(data_structure = data.frame(sp = tr$tip.label), phylogeny = list(sp = tr),
  #   n_pop = 300, seed = 1, parameters = p("sp"))
  # expect_equal(mean(sim_col(s, "g")), 5, tolerance = 0.05)

  K <- diag(5); K[1, 2] <- K[2, 1] <- 0.9; rownames(K) <- colnames(K) <- letters[1:5]
  s <- simulate_population(data_structure = data.frame(g = letters[1:5]), cov_str = list(g = K),
    n_pop = 2000, seed = 1, parameters = setNames(list(list(names = "u", mean = 5, vcov = 1, beta = 1), list(vcov = 1)), c("g", "residual")))
  expect_equal(unname(rowMeans(sim_col(s, "u"))), rep(5, 5), tolerance = 0.03)
})

test_that("cov_str can be larger than the data and is checked", {
  K <- diag(5); K[1, 2] <- K[2, 1] <- 0.9; rownames(K) <- colnames(K) <- letters[1:5]
  p <- list(g = list(names = "u", vcov = 1, beta = 1), residual = list(vcov = 1))
  s <- simulate_population(data_structure = data.frame(g = letters[1:3]), cov_str = list(g = K), parameters = p)
  expect_equal(nrow(get_population_data(s)), 3)
  K[1, 2] <- 0.1
  expect_error(simulate_population(data_structure = data.frame(g = letters[1:5]), cov_str = list(g = K), parameters = p), "symmetric")
})

test_that("covariates use the values in the data structure", {
  s <- simulate_population(data_structure = data.frame(time = rep(c(0, 5, 20), each = 10)), seed = 1,
    parameters = list(time = list(covariate = TRUE, beta = 1), residual = list(vcov = 1e-6)))
  d <- get_population_data(s)
  expect_equal(as.vector(tapply(d$y, d$time, mean)), c(0, 5, 20), tolerance = 0.01)
  expect_error(simulate_population(data_structure = data.frame(sex = rep(1:2, 5)),
    parameters = list(sex = list(fixed = TRUE, covariate = TRUE, beta = c(1, 2)), residual = list(vcov = 1))), "fixed")
})

test_that("get_parameters works with interactions", {
  s <- simulate_population(n = 10, parameters = list(observation = list(names = c("a", "b"), beta = c(1, 1)),
    interactions = list(names = "a:b", beta = 2), residual = list(vcov = 1)))
  expect_equal(unname(get_parameters(s)["a:b_beta"]), 2)
})

test_that("names on intercepts and betas are ignored (#32, #27)", {
  s <- simulate_population(n = 10, model = "y=intercept+residual", parameters = list(intercept = c(a = 1), residual = list(vcov = 1)))
  expect_s3_class(s, "squid")
  k <- data.frame(z = rnorm(10))
  s <- simulate_population(n = 10, n_response = 2, parameters = list(residual = list(vcov = c(1, 1))),
    known_predictors = list(predictors = k, beta = matrix(1:2, 1, 2, dimnames = list("z", c("p", "q")))))
  expect_equal(colnames(s$y[[1]]), c("y1", "y2"))
})

test_that("small and awkward data structures work", {
  expect_s3_class(simulate_population(n = 1, parameters = res), "squid")
  expect_s3_class(simulate_population(data_structure = make_structure("a(1)", repeat_obs = 5),
    parameters = list(a = list(vcov = 1), residual = list(vcov = 1))), "squid")
  s <- simulate_population(data_structure = data.frame(sex = c(1, 2, NA, 1, 2)),
    parameters = list(sex = list(fixed = TRUE, beta = c(1, 2)), residual = list(vcov = 1)))
  expect_equal(which(is.na(s$y[[1]])), 3)
})

test_that("invalid names are caught", {
  expect_error(simulate_population(n = 10, parameters = list(observation = list(names = "a^b", beta = 1), residual = list(vcov = 1))), "alphanumeric")
  expect_error(simulate_population(n = 10, parameters = list(observation = list(names = "y", beta = 1), residual = list(vcov = 1))), "unique")
  expect_error(simulate_population(n = 10, model = "y=a+residual", parameters = list(observation = list(names = "a", beta = 1, a_raw = 1), residual = list(vcov = 1))), "extra parameters")
})

test_that("make_structure gives informative errors", {
  expect_error(make_structure("individual"), "Could not read")
  expect_error(make_structure("a(0)"), "Could not read")
  expect_error(make_structure("a(2)+a(3)"), "unique")
  expect_error(make_structure("a(2)+a:b"), "Interactions")
})

test_that("nested sampling works with one factor and single level groups", {
  ds <- make_structure("individual(10)", repeat_obs = 4)
  s <- simulate_population(data_structure = ds, parameters = res, sample_type = "nested",
    sample_param = matrix(5, dimnames = list(NULL, "individual")))
  expect_equal(length(unique(get_sample_data(s)$individual)), 5)

  d <- data.frame(pop = rep(1:2, each = 8), ind = c(rep(9, 8), rep(1:4, each = 2)))
  n <- replicate(20, {
    s <- simulate_population(data_structure = d, parameters = res, sample_type = "nested",
      sample_param = matrix(c(2, 1), 1, 2, dimnames = list(NULL, c("pop", "ind"))))
    x <- get_sample_data(s)
    length(unique(x$ind[x$pop == 1]))
  })
  expect_true(all(n == 1))
})

test_that("missing sampling returns a list when samples are the same size", {
  s <- simulate_population(n = 20, parameters = res, sample_type = "missing", sample_param = c("100", "100"))
  expect_equal(nrow(get_sample_data(s, sample_set = 2)), 20)
})

test_that("temporal sampling covers all time points", {
  td <- make_structure("day(20)+individual(5)")
  s <- simulate_population(data_structure = td, parameters = res, sample_type = "temporal",
    sample_param = list(time = "day", group = "individual", variance = 0, n = 20))
  expect_equal(as.vector(table(get_sample_data(s)$individual)), rep(20, 5))
  r <- replicate(100, {
    s <- simulate_population(data_structure = td, parameters = res, sample_type = "temporal",
      sample_param = list(time = "day", group = "individual", variance = 0.5, n = 10))
    range(get_sample_data(s)$day)
  })
  expect_equal(range(r), c(1, 20))
})

test_that("get_sample_data complains if there are no samples", {
  expect_error(get_sample_data(simulate_population(n = 5, parameters = res)), "samples")
})

test_that("exp2lat and lat2exp are inverses and check input", {
  cv <- matrix(c(1, 0.2, 0.2, 2), 2)
  r <- exp2lat(c(2, 3), cv)
  b <- lat2exp(r$mean, r$cov)
  expect_equal(b$mean, c(2, 3))
  expect_equal(b$cov, cv)
  expect_error(exp2lat(1, "a"), "numeric")
  expect_error(exp2lat(1, -1), ">=0")
})

test_that("#33 can be done with an interaction", {
  d <- make_structure("population(2)/individual(3)+year(2)+population:year", int_names = "pop_year")
  expect_equal(length(unique(d$pop_year)), 4)
  expect_equal(nrow(unique(d[, c("population", "pop_year")])), 4)
})


test_that("pedigree handles singular G and extra columns", {
  s <- simulate_population(data_structure = data.frame(animal = ped$id), pedigree = list(animal = cbind(ped, sex = 1)),
    parameters = list(animal = list(names = c("g1", "g2"), vcov = matrix(1, 2, 2), beta = c(1, 1)), residual = list(vcov = 1)))
  d <- get_population_data(s)
  expect_equal(d$g1, d$g2)
})

test_that("interactions work with numeric names", {
  s <- simulate_population(data_structure = data.frame(yr = rep(c(2001, 2002), each = 5)), parameters = list(
    yr = list(fixed = TRUE, names = c("2001", "2002"), beta = c(0, 0)),
    observation = list(names = "x", beta = 0), interactions = list(names = "2001:x", beta = 1), residual = list(vcov = 1)))
  d <- get_population_data(s)
  expect_equal(d[["2001:x"]], d$x * (d$yr == 2001))
})


test_that("temporal sampling handles non-consecutive and named times", {
  d <- make_structure("day(10)+individual(3)"); d$day <- d$day * 2
  s <- simulate_population(data_structure = d, parameters = res, sample_type = "temporal",
    sample_param = list(time = "day", group = "individual", variance = 0, n = 10))
  expect_equal(nrow(get_sample_data(s)), 30)
  d <- make_structure("day(4)+individual(3)", level_names = list(day = letters[1:4]))
  s <- simulate_population(data_structure = d, parameters = res, sample_type = "temporal",
    sample_param = list(time = "day", group = "individual", variance = 0, n = 4))
  expect_equal(nrow(get_sample_data(s)), 12)
})

test_that("nested sampling handles factors, no data structure and bad requests", {
  d <- make_structure("pop(3)/ind(4)", repeat_obs = 2)
  d$pop <- factor(paste0("p", d$pop)); d$ind <- factor(paste0("i", d$ind))
  p <- matrix(c(2, 1), 1, 2, dimnames = list(NULL, c("pop", "ind")))
  s <- simulate_population(data_structure = d, parameters = res, sample_type = "nested", sample_param = p)
  expect_equal(nrow(get_sample_data(s)), 4)
  s <- simulate_population(n = 10, parameters = res, sample_type = "nested", sample_param = matrix(2, dimnames = list(NULL, "observation")))
  expect_equal(nrow(get_sample_data(s)), 2)
  expect_error(simulate_population(data_structure = d, parameters = res, sample_type = "nested", sample_param = matrix(2)), "column names")
  p[1, 2] <- 9
  expect_error(simulate_population(data_structure = d, parameters = res, sample_type = "nested", sample_param = p), "maximum available")
})

test_that("model must assign a response", {
  expect_error(simulate_population(n = 5, model = "a+residual", parameters = list(
    observation = list(names = "a", beta = 1), residual = list(vcov = 1))), "assign")
})

test_that("nested sampling warns when some groups are short", {
  d <- data.frame(pop = rep(1:2, c(2, 8)), ind = c(1, 1, rep(2:5, each = 2)))
  expect_warning(simulate_population(data_structure = d, parameters = res, sample_type = "nested",
    sample_param = matrix(c(2, 2), 1, 2, dimnames = list(NULL, c("pop", "ind")))), "fewer than")
})
