
generate_internal_structure <- function(data_structure, n, parameters, n_response=1, response_names, known_predictors, model, index_link, family="gaussian", link="identity", pedigree, phylogeny, cov_str,sample_type, sample_param, sample_plot=FALSE, n_pop=1, seed, verbose=FALSE,suppress_index_warning=FALSE){
  
  if(missing("n") & missing("data_structure")){
    stop("Either 'n' or 'data_structure' need to be specified")
  }else if(missing("n")){
    n <- nrow(data_structure)
  }else if(!missing("n") & !missing("data_structure")){
    if(nrow(data_structure)!=n) stop("'n' and nrow(data_structure) are not equal. Only one needs to be specified.")
  }
  
  if(missing("parameters")) stop("'parameters' need to be specified")

  if(!all(link %in% c("identity", "log", "inverse", "sqrt", "logit", "probit", "cloglog"))) stop("Link must be 'identity', 'log', 'inverse', 'sqrt', 'logit', 'probit', 'cloglog'")
  if(!(length(link)==n_response || length(link)==1))  stop("Link must either be length 1 or same length as the number of responses")
  
  if(!all(family %in% c("gaussian", "poisson", "binomial"))) stop("Family must be 'gaussian', 'poisson', 'binomial'")
  if(!(length(family)==n_response || length(family)==1))  stop("Family must either be length 1 or same length as the number of responses")


  if(n_response > 1 & !missing("model")) stop("Currently cannot specify multiple responses and a model formula")

  if(!missing(response_names)) {
    if(!missing("model")) message("response_names is ignored when a model formula is specified") 
    if(length(response_names)!=n_response) stop("response_names needs to be the same length as n_response")  
  }

  # if(!missing(pedigree)){
  #   if(missing(pedigree_type)){
  #     pedigree_type <- as.list(rep("A",length(pedigree)))
  #     names(pedigree_type) <- names(pedigree)
  #   }else{
  #     if(!pedigree_type %in% c("A","D","E"))stop("pedigree_type must be either 'A','D' or 'E'")
  #     if(length(pedigree)!=length(pedigree_type)) stop("pedigree and pedigree_type need to be the same length")
  #     if(sort(names(pedigree))==sort(names(pedigree_type))) stop("names of pedigree and pedigree_type need to match")
  #   }
  # }
  # if(!missing(phylogeny)){
  #   if(missing(phylogeny_type)){
  #     phylogeny_type <- as.list(rep("brownian",length(phylogeny)))
  #     names(phylogeny_type) <- names(phylogeny)
  #   }else{
  #     if(!phylogeny_type %in% c("brownian","OU")) stop("phylogeny_type must be either 'brownian' or 'OU'")
  #     if(length(phylogeny)!=length(phylogeny_type)) stop("phylogeny and phylogeny_type need to be the same length")
  #     if(sort(names(phylogeny))==sort(names(phylogeny_type))) stop("names of phylogeny and phylogeny_type need to match")
  #   }
  # }

  ## gets the arguments into a list that is added to for the output
  lapply(as.list(environment()), function(x) if (!is.list(x) &&length(x)==1 && x=="") NULL else x)
}

# simulate_population <- function(data_structure, n, parameters, n_response=1, response_names, known_predictors, model, index_link, family="gaussian", link="identity", pedigree, pedigree_type, phylogeny, phylogeny_type, cov_str,sample_type, sample_param, sample_plot=FALSE, n_pop=1, verbose=FALSE){


#   # input<- lapply(as.list(environment()), function(x) if (!is.list(x) &&length(x)==1 && x=="") NULL else x)
#   do.call(generate_internal_structure,as.list(environment()))

# }

