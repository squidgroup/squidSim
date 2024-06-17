
#' @title simulate_population
#' @description Simulate population level data
#' @param data_structure A matrix or data.frame with a named column for each grouping factor, including the levels
#' @param n Sample size when data_structure is not specified
#' @param parameters A list of parameters for each hierarchical level. See details.
#' @param n_response The number of response variables, defaults to 1. 
#' @param response_names Names given to response variables. Defaults to 'y', or c('y1','y2',...) if n_response>1. Not used if model is specified. 
#' @param known_predictors This argument provides a way of inputting existing predictor variables. This argument takes a list, with items 'predictors' and 'betas', where 'predictors' is a matrix or data.frame, the number of rows of which MUST equal 'n' or the number of rows in 'data_structure'. 
#' @param model Optional. A formula description of the simulation model. See details.
#' @param index_link Optional. Make new factors in the data structure indexed by other factors. This is only really used in the context of the model argument. Takes a list, the names of which indicate the new grouping factor name and the elements of the list represent which grouping factors should be linked, for example "mother-ID" would index mother with the same indexes as ID.  
#' @param family A description of the error distribution. Options are 'gaussian' (default), 'poisson' and 'binomial'. 'binomial' generates a binary response variable.
#' @param link A description of the link function distribution. Options are 'identity' (default),'log', 'inverse', 'sqrt', 'logit' and 'probit'.
#' @param pedigree A list of pedigrees for each hierarchical level. Each pedigree must be matrix or data.frame, that is at least 3 columns, which correspond to ID, dam and sire. The name in the pedigree list must match a name in the parameter list.
#' @param pedigree_type A list describing what kind of genetic variance is to be simulated from each pedigree. Default is 'A', other options are 'D' (dominance) and 'E' (epistatic). Makes use of relationship matrices created by the MCMCglmm and nadiv packages.
#' @param phylogeny A list of phylogenies for each hierarchical level. Each phylogeny should be phylo class. The name in the phylogeny list must match a name in the parameter list.
#' @param phylogeny_type A list describing what mode of evolution should be simulated from each phylogeny. Options are 'brownian'(default) or 'OU'. 
#' @param cov_str A list of covariance structures for each hierarchical level. The name in the cov_str list must match a name in the parameter list.
#' @param sample_type Type of sampling, must be one of 'nested', 'missing' or temporal. If not specified, then no sampling is done. See details
#' @param sample_param A set of parameters, specific to the sampling type. See details.
#' @param sample_plot Logical. Should illustrative plots be made - defaults to FALSE - currently not implemented.
#' @param n_pop Number of populations. Default = 1
#' @param verbose Logical. Whether to print diagnostics. Useful for debugging. Defaults to FALSE
#' @details 
#' A detailed vignette can be found at http://squidgroup.org/squidSim_vignette/
#' 
#' The parameters list contains one (or more) list for each hierarchical level that you want to simulate at. A residual list is always need, specifying variances/covariances for the residual. Additionally, the parameter list can also be provided with an intercept vector and interactions list. For each item in the parameter list (excluding intercept, interactions, and residual), the following can be specified (but all have default values): 
#' names - vector containing the names of predictors simulated at this level
#'  group - character string relating to the data_structure
#' mean - vector of means for the predictor variables
#' vcov or vcorr - Either a vector of variances, or a variance-covariance/correlation matrix, for the predictor variables
#' beta - vector of effect sizes (or matrix with n_response columns when n_response>1)
#' fixed - Logical, indicating whether the effects for the levels are fixed or to be simulated
#' covariate - Logical, indicating whether the indexes in the data structure are to be used as a continuous variable rather than simulating one
#' functions - vector - transformation to be applied to the response variable. Defaults to ‘identity’.
#' A more detailed explanation can be found at http://squidgroup.org/squidSim_vignette/1.9-parameter-list-summary.html 
#' 
#' The model argument is character string which explicitly tells the simulate_population function how to put the simulated predictors together to form the response variable. For example, if the predictors temperature and rainfall had been specified in the parameter list, providing the model argument with "y = temperature + rainfall + residual" would result in generating a response variable 'y' in the same way the siumulate_population function does by default. For more detailed information see http://squidgroup.org/squidSim_vignette/1.7-modeleq.html
#' 
#' Different sampling schemes can be implemented, (sample_type can be 'nested', 'missing' or temporal). The sample_param takes a different form depending on the sample_type. See http://squidgroup.org/squidSim_vignette/7-sampling.html for full details.
#' 
#' 
#' @author Joel Pick - joel.l.pick@gmail.com
#' @return a squid object, which is a list including all inputs and simulated data.
#' @examples
#' # simple linear model with three predictors variables
#' squid_data <- simulate_population(
#'   n=50,
#'   parameters = list(
#'     observation = list(
#'       names = c("temperature","rainfall", "wind"),
#'       beta = c(0.5,-0.3, 0.4)    
#'   ),
#'     residual = list(
#'       vcov = 1
#'     )
#'   )
#' )
#' 
#' @export
simulate_population <- function(data_structure, n, parameters, n_response=1, response_names, known_predictors, model, index_link, family="gaussian", link="identity", pedigree, pedigree_type, phylogeny, phylogeny_type, cov_str,sample_type, sample_param, sample_plot=FALSE, n_pop=1, verbose=FALSE){

  if(verbose) cat("checking input\n")

  if(!all(link %in% c("identity", "log", "inverse", "sqrt", "logit", "probit"))) stop("Link must be 'identity', 'log', 'inverse', 'sqrt', 'logit', 'probit'")
  if(!(length(link)==n_response || length(link)==1))  stop("Link must either be length 1 or same length as the number of responses")
  
  if(!all(family %in% c("gaussian", "poisson", "binomial"))) stop("Family must be 'gaussian', 'poisson', 'binomial'")
  if(!(length(family)==n_response || length(family)==1))  stop("Family must either be length 1 or same length as the number of responses")

  if(missing("n") & missing("data_structure")){
    stop("Either 'n' or 'data_structure' need to be specified")
  }else if(missing("n")){
    n <- nrow(data_structure)
  }else if(!missing("n") & !missing("data_structure")){
    if(nrow(data_structure)!=n) stop("'n' and nrow(data_structure) are not equal. Only one needs to be specified.")
  }
  
  if(!missing("known_predictors")){
    if(n!=nrow(known_predictors[["predictors"]])) stop("The number of observation specified does not match the number of rows in known_predictors")
  }

  if(n_response > 1 & !missing("model")) stop("Currently cannot specify multiple responses and a model formula")

  if(!missing(response_names)) {
    if(!missing("model")) message("response_names is ignored when a model formula is specified") 
    if(length(response_names)!=n_response) stop("response_names needs to be the same length as n_response")  
  }

  if(!missing(pedigree)){
    if(missing(pedigree_type)){
      pedigree_type <- as.list(rep("A",length(pedigree)))
      names(pedigree_type) <- names(pedigree)
    }else{
      if(!pedigree_type %in% c("A","D","E"))stop("pedigree_type must be either 'A','D' or 'E'")
      if(length(pedigree)!=length(pedigree_type)) stop("pedigree and pedigree_type need to be the same length")
      if(sort(names(pedigree))==sort(names(pedigree))) stop("names of pedigree and pedigree_type need to match")
    }
  }
  if(!missing(phylogeny)){
    if(missing(phylogeny_type)){
      phylogeny_type <- as.list(rep("brownian",length(phylogeny)))
      names(phylogeny_type) <- names(phylogeny)
    }else{
      if(!phylogeny_type %in% c("brownian","OU")) stop("phylogeny_type must be either 'brownian' or 'OU'")
      if(length(phylogeny)!=length(phylogeny_type)) stop("phylogeny and phylogeny_type need to be the same length")
      if(sort(names(phylogeny))==sort(names(phylogeny))) stop("names of phylogeny and phylogeny_type need to match")
    }
  }






  ## gets the arguments into a list that is added to for the output
  output <- lapply(as.list(environment()), function(x) if (!is.list(x) &&length(x)==1 && x=="") NULL else x)

#####################  
###---Fill in parameter lists 
##################### 


  if(!missing(known_predictors)) {
    if(verbose) cat("checking known_predictors\n")
    output$known_predictors <- do.call(fill_preds, output)
  }

  if(verbose) cat("checking parameter list\n")

  output$parameters <- do.call(fill_parameters, output)

  # if(verbose) print(output$parameters)
  


#####################  
###---data structure as indexes
      ## index data_structure

#####################  
  
  if(verbose) cat("indexing data_structure\n")
  output$str_index <- do.call(index_factors, output)
  
  # if(verbose)  {
  #   print(head(output$str_index))
  #   print(tail(output$str_index))
  # }
#####################  
###---cov structures
#####################  


  ## check pedigree levels match data structure levels
  ## make function - that can check ped,phylo and covs
  
  # if(verbose) cat("checking pedigree/phylogeny \n")



  if(verbose) cat("generating covariance structures \n")

  output$cov_str_all <- do.call(cov_str_list, output)
  ## make cov_str with everything, then return it back to cov_str after predictors



#####################  
###---PREDICTORS 
#####################  
  if(verbose) cat("simulating predictors\n")

  output$predictors <- lapply(1:n_pop, function(x) do.call(sim_predictors, output))
  # output$predictors <- lapply(1:n_pop, function(x) cbind(do.call(sim_predictors, output), known_predictors))
  ## returns list of predictor matrices

  # output$cov_str <- cov_str


#####################  
###---GENERATE Y
##################### 
  if(verbose) cat("generating y\n")

  # y <- do.call(generate_y, output)
  y <- do.call(generate_y_list, output)


#####################  
###---TRANSFORM Y 
##################### 
  if(verbose) cat("transforming y\n")

  output$y <- lapply(y, function(x) transform_dist(x, family, link, output$response_names))


#####################  
###--- DO SAMPLING 
##################### 

  if(verbose) cat("sampling\n")

  if(!is.null(output$sample_type)) output$samples <- sample_population(output)


### work out what to return 

  output <- output[!names(output)%in%c("cov_str_all","str_index")]


  class(output) <- 'squid'
  return(output)
}

## problem that by default the predictors and the level IDs will have the same names
## - maybe append "_effects"


#' @title print.squid
#' @description Print method for class 'squid'
#' @param x an R object of class 'squid'
#' @param ... further arguments passed to or from other methods.
#' @export
print.squid <- function(x, ...){
  cat("Data simulated using squid \n
              /\\             
             /  \\            
            / /\\ \\           
            \\/  \\/            
            /    \\           
            |    |          
            |    |          
     0      |    |      0     
     /      \\____/      \\    
    {     __/(  )\\__     }   
     \\___/__\\_\\/_/__\\___/    
      / / / /    \\ \\ \\ \\     
     / / / {      } \\ \\ \\    
    { { /   \\    /   \\ } }   
    }  \\     0  0     /  {   
 0_/  { \\_0        0_/ }  \\_0
       \\              /      
        }            {       
       /              \\      
      0                0      
")
}



#' @title summary.squid
#' @description summary method for class 'squid'
#' @param object an R object of class 'squid'
#' @param ... further arguments passed to or from other methods.
#' @export
summary.squid <- function(object, ...){
  
  ## description of sampling.
}


#' @title get_population_data
#' @description Extracts population level data from a squid object
#' @param x an R object of class 'squid'
#' @param list Logical - whether to return data as a list or data_table (FALSE; default).
#' @param ... further arguments passed to or from other methods.
#' @export
get_population_data <- function(x,list=FALSE,...){

  # pop_list <- lapply(1:x$n_pop,function(i) data.table::data.table(cbind(x$y[[i]],x$predictors[[i]],x$data_structure,squid_pop=i)))
  pop_list <- lapply(1:x$n_pop,function(i) as.data.frame(cbind(x$y[[i]],x$predictors[[i]],x$data_structure,squid_pop=i)))
  if(list){
    return(pop_list)
  }else{
    do.call(rbind,pop_list)
  }

}


#' @title get_parameters
#' @description Extracts population level data from a squid object
#' @param x an R object of class 'squid'
#' @export

get_parameters <- function(x){
  intercept <- x$parameters$intercept
  param <- x$parameters[names(x$parameters)!="intercept"]
  p_out<- lapply(param,function(y){

    means <- y$mean
    names(means) <- paste0(y$names,"_mean")

    vars <- diag(y$vcov)
    names(vars) <- paste0(y$names,"_var")

    if(ncol(y$vcov)>1){ 
      covs <- y$vcov[lower.tri(y$vcov)]
      name_ind<-which(lower.tri(y$vcov), arr.ind=TRUE)
      names(covs) <- paste0(y$names[name_ind[,2]],":",y$names[name_ind[,1]],"_cov")
    }else{
      covs<-NULL
    }

    betas <- as.vector(y$beta)
    names(betas) <- if(x$n_response==1){ 
      paste0(y$names,"_beta") 
    }else{ 
      paste0(rep(y$names,x$n_response),"_beta_y",rep(1:x$n_response, each=nrow(y$beta) )) 
    }

    c(means,vars,covs,betas)  
  })
  names(p_out)<-NULL
  c(intercept,p_out, recursive=TRUE)

  ## extra_parameters
}