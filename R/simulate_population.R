
#' @title simulate_population
#' @description Simulate population level data
#' @param data_structure A matrix or dataframe with a named column for each grouping factor, including the levels
#' @param N Sample size when data_structure is not specified
#' @param parameters A list of parameters for each hierarchical level. See details.
#' @param N_response The number of response variables, defaults to 1. 
#' @param response_names Names given to response variables, defaults to y, or y1/y2... if there are multiple responses. Not used if model is specified. 
#' @param known_predictors This option provides a way of inputting existing predictor variables, without simulating all predictors. This argument takes a list, with item 'predictors' and 'betas', where 'predictors' is a matrix or dataframe 
#' @param model Optional. 
#' @param family A description of the error distribution. Default "gaussian".
#' @param link A description of the link function distribution. Default "identity".
#' @param pedigree A list of pedigrees for each hierarchical level. Each pedigree must be matrix or data.frame, that is at least 3 columns, which correspond to ID, dam and sire.
#' @param pedigree_type A list describing what kind of genetic variance is to be simulated from each pedigree. Default is 'additive', other options are 'dominance' and 'epistatic'. Makes use of relationship matrices created by the nadiv package.
#' @param phylogeny A list of phylogenies for each hierarchical level. Each pedigree should be phylo class.
#' @param phylogeny_type A list describing what mode of evolution should be simulated from each phylogeny. Options are 'brownian'(default) or 'OU'. 
#' @param cov_str A list of covariance structures for each hierarchical level. 
#' @param sample_type Type of sampling, needs to be one of 'nested', 'missing' or temporal. See details
#' @param sample_param A set of parameters, specific to the sampling type. See details.
#' @param sample_plot Logical. Should illustrative plots be made - defaults to FALSE.
#' @param N_pop Number of populations. Default = 1
#' @param verbose Logical. Whether to print diagnostics - defaults to FALSE
#' @details Parameter list ... 
#' @author Joel Pick - joel.l.pick@gmail.com
#' @return a squid object, which is a list including all inputs and simulated data.
#' @examples
#' # simple linear model with three predictors variables
#' squid_data <- simulate_population(
#'   parameters = list(
#'     observation = list(
#'       names = c("temperature","rainfall", "wind"),
#'       beta = c(0.5,-0.3, 0.4)    
#'   ),
#'     residual = list(
#'       vcov = 1
#'     )
#'   ),
#'   N=2000
#' )
#' 
#' @export
simulate_population <- function(data_structure, N, parameters, N_response=1, response_names, known_predictors, model, family="gaussian", link="identity", pedigree, pedigree_type, phylogeny, phylogeny_type, cov_str,sample_type, sample_param, sample_plot=FALSE, N_pop=1, verbose=FALSE){

  if(verbose) cat("checking input\n")

  if(!all(link %in% c("identity", "log", "inverse", "sqrt", "logit", "probit"))) stop("Link must be 'identity', 'log', 'inverse', 'sqrt', 'logit', 'probit'")
  if(!(length(link)==N_response || length(link)==1))  stop("Link must either be length 1 or same length as the number of responses")
  
  if(!all(family %in% c("gaussian", "poisson", "binomial"))) stop("Family must be 'gaussian', 'poisson', 'binomial'")
  if(!(length(family)==N_response || length(family)==1))  stop("Family must either be length 1 or same length as the number of responses")

  if(missing("N") & missing("data_structure")){
    stop("Either N or data_structure need to be specified")
  }else if(missing("N")){
    N <- nrow(data_structure)
  }else if(!missing("N") & !missing("data_structure")){
    if(nrow(data_structure)!=N) stop("N and nrow(data_structure) are not equal. Only one needs to be specified.")
  }
  
  if(!missing("known_predictors")){
    if(N!=nrow(known_predictors[["predictors"]])) stop("The number of observation specified does not match the number of rows in known_predictors")
  }

  if(N_response > 1 & !missing("model")) stop("Currently cannot specify multiple responses and a model formula")

  if(!missing(response_names)) {
    if(!missing("model")) message("response_names is ignored when a model formula is specified") 
    if(length(response_names)!=N_response) stop("response_names needs to be the same length as N_response")  
  }


  ## gets the arguments into a list that is added to for the output
  output <- lapply(as.list(environment()), function(x) if (length(x)==1 && x=="") NULL else x)

#####################  
###---Fill in parameter lists 
##################### 


  if(verbose) cat("checking parameter list\n")

  output$parameters <- do.call(fill_parameters, output)


  if(verbose) cat("checking known_predictors\n")

  if(!missing(known_predictors)) output$known_predictors <- do.call(fill_preds, output)
  


#####################  
###---data structure as indexes
      ## index data_structure

#####################  
  
  if(verbose) cat("indexing data_structure\n")
  output$str_index <- do.call(index_factors, output)
  

#####################  
###---cov structures
#####################  


  ## check pedigree levels match data structure levels
  ## make function - that can check ped,phylo and covs
  
  if(verbose) cat("checking pedigree/phylogeny \n")

  if(!missing(pedigree)){
    if(missing(pedigree_type)){
      output$pedigree_type <- as.list(rep("additive",length(pedigree)))
      names(output$pedigree_type) <- names(pedigree)
    }else{
      if(!pedigree_type %in% c("additive","dominance","epistatic"))stop("phylogeny_type must be wither 'additive','dominance' or 'epistatic'")
      if(length(pedigree)!=length(pedigree_type)) stop("pedigree and pedigree_type need to be the same length")
      if(sort(names(pedigree))==sort(names(pedigree))) stop("names of pedigree and pedigree_type need to match")
    }
  }
  if(!missing(phylogeny)){
    if(missing(phylogeny_type)){
      output$phylogeny_type <- as.list(rep("brownian",length(phylogeny)))
      names(output$phylogeny_type) <- names(phylogeny)
    }else{
      if(!phylogeny_type %in% c("brownian","OU")) stop("phylogeny_type must be wither 'brownian' or 'OU'")
      if(length(phylogeny)!=length(phylogeny_type)) stop("phylogeny and phylogeny_type need to be the same length")
      if(sort(names(phylogeny))==sort(names(phylogeny))) stop("names of phylogeny and phylogeny_type need to match")
    }
  }

  if(verbose) cat("generating covariance structures \n")

  output$cov_str_all <- do.call(cov_str_list, output)
  ## make cov_str with everything, then return it back to cov_str after predictors



#####################  
###---PREDICTORS 
#####################  
  if(verbose) cat("simulating predictors\n")

  output$predictors <- lapply(1:N_pop, function(x) do.call(sim_predictors, output))
  # output$predictors <- lapply(1:N_pop, function(x) cbind(do.call(sim_predictors, output), known_predictors))
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

  # data.table(cbind(x$y,x$predictors,x$data_structure))

  pop_list <- lapply(1:x$N_pop,function(i) data.table::data.table(cbind(x$y[[i]],x$predictors[[i]],x$data_structure,squid_pop=i)))

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
    names(betas) <- if(x$N_response==1){ 
      paste0(y$names,"_beta") 
    }else{ 
      paste0(rep(y$names,x$N_response),"_beta_y",rep(1:x$N_response, each=nrow(y$beta) )) 
    }

    c(means,vars,covs,betas)  
  })
  names(p_out)<-NULL
  c(intercept,p_out, recursive=TRUE)

  ## extra_parameters
}