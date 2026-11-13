
#' @title simulate_population
#' @description Simulate population level data
#' @param data_structure A matrix or data.frame with a named column for each grouping factor, including the levels
#' @param n Sample size when data_structure is not specified
#' @param parameters A list of parameters for each hierarchical level. See details.
#' @param n_response The number of response variables, defaults to 1. 
#' @param response_names Names given to response variables. Defaults to 'y', or c('y1','y2',...) if n_response>1. Not used if model is specified. 
#' @param known_predictors This argument provides a way of inputting existing predictor variables. This argument takes a list, with items 'predictors' and 'beta', where 'predictors' is a matrix or data.frame, the number of rows of which MUST equal 'n' or the number of rows in 'data_structure'. The column names of 'predictors' are used as the variable names, and the 'beta's are assumed to be in the same order as the predictors.
#' @param model Optional. A equation description of the simulation model given as a character string. See details.
#' @param index_link Optional. Make new factors in the data structure indexed by other factors. This is only really used in the context of the model argument. Takes a list, the names of which indicate the new grouping factor name and the elements of the list represent which grouping factors should be linked, for example "mother-ID" would index mother with the same indexes as ID.
#' @param family A description of the error distribution. Options are 'gaussian' (default), 'poisson' and 'binomial'. 'binomial' generates a binary response variable.
#' @param link A description of the link function distribution. Options are 'identity' (default),'log', 'inverse', 'sqrt', 'logit', 'probit' and 'cloglog'.
#' @param pedigree A list of pedigrees for each hierarchical level. Each pedigree must be matrix or data.frame, that is at least 3 columns, which correspond to ID, dam and sire. The name in the pedigree list must match a name in the parameter list.
#' @param phylogeny A list of phylogenies for each hierarchical level. Each phylogeny should be phylo class. The name in the phylogeny list must match a name in the parameter list.
#' @param cov_str A list of covariance structures for each hierarchical level. The name in the cov_str list must match a name in the parameter list.
#' @param sample_type Type of sampling, must be one of 'nested', 'missing', 'survival' or 'temporal.' If not specified, then no sampling is done. See details
#' @param sample_param A set of parameters, specific to the sampling type. See details.
#' @param n_pop Number of populations. Default = 1
#' @param seed Seed for random number generation. Setting seed allows for simulations to be easily reproduced. If seed is not set, one is set internally, which is then saved and returned in the squidSim object
#' @param verbose Logical. Whether to print diagnostics. Useful for debugging. Defaults to FALSE
#' @param suppress_index_warning Logical. Whether to print warnings relating to the index-link argument. Useful to switch off if using in large number of simulations. Defaults to FALSE.
#' @details 
#' A detailed vignette can be found at https://squidgroup.org/squidSim_vignette/
#' 
#' The parameters list contains one (or more) list for each hierarchical level that you want to simulate at. A residual list is always needed, specifying variances/covariances for the residual. Additionally, the parameter list can also be provided with an intercept vector and interactions list. For each item in the parameter list (excluding intercept, interactions, and residual), the following can be specified (but all have default values): 
#' names - vector containing the names of predictors simulated at this level
#'  group - character string relating to the data_structure
#' mean - vector of means for the predictor variables
#' vcov or vcorr - Either a vector of variances, or a variance-covariance/correlation matrix, for the predictor variables
#' beta - vector of effect sizes (or matrix with n_response columns when n_response>1)
#' fixed - Logical, indicating whether the effects for the levels are fixed or to be simulated
#' covariate - Logical, indicating whether the indexes in the data structure are to be used as a continuous variable rather than simulating one
#' functions - vector - transformation to be applied to the response variable. Defaults to ‘identity’.
#' A more detailed explanation can be found at https://squidgroup.org/squidSim_vignette/1.9-parameter-list-summary.html 
#' 
#' The model argument is character string which explicitly tells the simulate_population function how to put the simulated predictors together to form the response variable. For example, if the predictors temperature and rainfall had been specified in the parameter list, providing the model argument with "y = temperature + rainfall + residual" would result in generating a response variable 'y' in the same way the simulate_population function does by default. For more detailed information see http://squidgroup.org/squidSim_vignette/1.7-modeleq.html
#' 
#' Different sampling schemes can be implemented, (sample_type can be 'nested', 'missing', 'survival' or 'temporal'). The sample_param takes a different form depending on the sample_type. See https://squidgroup.org/squidSim_vignette/7-sampling.html for full details.
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
simulate_population <- function(data_structure, n, parameters, n_response=1, response_names, known_predictors, model, index_link, family="gaussian", link="identity", pedigree, phylogeny, cov_str,sample_type, sample_param, n_pop=1, seed, verbose=FALSE,suppress_index_warning=FALSE){


  if(verbose) cat("checking input\n")
  output <- do.call(generate_internal_structure,as.list(environment()))
  

#####################  
###---Fill in parameter lists 
##################### 


  # if(!missing(known_predictors)) {
  if(!is.null(output$known_predictors)) {
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

  # if(is.null(output$seed))output$seed<-sample(1:10000,1)
## if no seed is given, draw one from the user's RNG stream. The stream is restored on exit, so that repeated calls keep moving along it, rather than being reset to one of a small set of seeds, and set.seed() here has no side effects on the user's session
 
  # if not RNG steam, initiates RNG steam
  if(!exists(".Random.seed", envir=globalenv(), inherits=FALSE)) stats::runif(1)
  
  # if no seed is given, draw one from the user's RNG stream.
  if(is.null(output$seed)) output$seed <- sample.int(.Machine$integer.max,1)
  
  # get RNG steam
  old_seed <- get(".Random.seed", envir=globalenv(), inherits=FALSE)

  # restore RNG steam on exit, so that repeated calls keep moving along it, rather than being reset to one of a small set of seeds, and set.seed() here has no side effects on the user's session
  on.exit(assign(".Random.seed", old_seed, envir=globalenv()), add=TRUE)

  set.seed(output$seed)
  
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

  if(!is.null(output$sample_type) & !is.null(output$sample_param)) output$samples <- sample_population(output)

  if((is.null(output$sample_type) & !is.null(output$sample_param)) | (!is.null(output$sample_type) & is.null(output$sample_param))) warning("sample_type and sample_param both need to be specified for sampling to be done. Simulations have been returned without sampling", call.=FALSE)

### work out what to return 

  output <- output[!names(output)%in%c("cov_str_all","str_index")]


  class(output) <- 'squid'
  return(output)
}



#' @title print.squid
#' @description Print method for class 'squid'
#' @param x an R object of class 'squid'
#' @param ... further arguments passed to or from other methods.
#' @return x, invisibly
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
  invisible(x)
}


#' @title summary.squid
#' @description summary method for class 'squid'
#' @param object an R object of class 'squid'
#' @param ... further arguments passed to or from other methods.
#' @return A list, printed as a side effect, with the number of populations, observations and responses, and the sampling type.
#' @export
summary.squid <- function(object, ...){
  
  out <- list(n_pop=object$n_pop, n=object$n, n_response=object$n_response, seed=object$seed, sample_type=object$sample_type)
  cat("squid simulation:",out$n_pop,"population(s),",out$n,"observations,",out$n_response,"response(s)\n")
  cat("Sampling:",if(is.null(out$sample_type)) "none" else out$sample_type,"\n")
  invisible(out)
}


#' @title get_population_data
#' @description Extracts population level data from a squid object
#' @param x an R object of class 'squid'
#' @param list Logical - whether to return data as a list or data.frame (FALSE; default).
#' @param ... further arguments passed to or from other methods.
#' @return A data.frame with the response(s), predictors and data structure for all populations, or a list of data.frames, one for each population, if list=TRUE.
#' @examples
#' squid_data <- simulate_population(n=10, parameters=list(residual=list(vcov=1)))
#' get_population_data(squid_data)
#' @export
get_population_data <- function(x,list=FALSE,...){

  pop_list <- lapply(1:x$n_pop,function(i) as.data.frame(cbind(x$y[[i]],x$predictors[[i]],as.data.frame(x$data_structure),squid_pop=i)))
  if(list){
    return(pop_list)
  }else{
    do.call(rbind,pop_list)
  }

}


#' @title get_parameters
#' @description Extracts the simulation parameters from a squid object
#' @param x an R object of class 'squid'
#' @return A named vector of the means, variances, covariances and betas used in the simulation.
#' @examples
#' squid_data <- simulate_population(n=10, parameters=list(residual=list(vcov=1)))
#' get_parameters(squid_data)
#' @export

get_parameters <- function(x){
  intercept <- x$parameters$intercept
  param <- x$parameters[names(x$parameters)!="intercept"]
  names(intercept) <- if(length(intercept)==1) "intercept" else paste0("intercept_y",seq_along(intercept))

  p_out<- lapply(param,function(y){
  ## interactions only have names and betas
    if(is.null(y$vcov)){
      means <- vars <- covs <- NULL
    }else{
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
  known <- x$known_predictors$beta
  known_betas <- as.vector(known)
  if(!is.null(known)) names(known_betas) <- if(x$n_response==1){ 
    paste0(rownames(known),"_beta") 
  }else{ 
    paste0(rep(rownames(known),x$n_response),"_beta_y",rep(1:x$n_response, each=nrow(known) )) 
  }

  c(intercept,p_out,known_betas, recursive=TRUE)

  ## extra_parameters
}