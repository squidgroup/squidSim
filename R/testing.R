# devtools::load_all("./github/squidSim/R")

# ## simplest model
# output<- generate_internal_structure(n=5,parameters=list(residual=list(vcov=1)))

# if(!is.null(output$known_predictors)) {   
#   output$known_predictors <- do.call(fill_preds, output)
# }

# output$parameters <- do.call(fill_parameters, output)

# output$str_index <- do.call(index_factors, output)
#    