
n_phenotypes <- function(parameters){
  j <- sapply(parameters[names(parameters) !="interactions"], function(x) ncol(x[["beta"]]))
  if(length(unique(j))!= 1) stop("The number of phenotypes (columns in beta) are not consistent across hierarchical levels in the parameter list", call.=FALSE)
  return(unique(j))
}
# [names(parameters) !="interactions"]
# n_phenotypes(parameters)



## I've used loops rather than apply functions in here because then the original parameter list can then be added to rather than new lists made - this will be slightly slower but very negligible given their size
fill_parameters <- function(parameters,data_structure, n, n_response, response_names,known_predictors,...){

  reserved_p_names <- c("intercept","observation","residual","interactions") 

  # Check whether list is given
  if(!is.list(parameters)) stop("parameters are not provided as a list", call.=FALSE)

  # 

#########
## Group Names
#########
  
  # Check whether group specified -If not, then give name in list
  # Check whether length(group)==1 - If not, give warning and use only first group
 
  #i <- names(parameters)[1]
  for (i in names(parameters)[!names(parameters) %in% c("intercept","interactions")]){
    group <-parameters[[i]][["group"]]
    if(is.null(group)) group <- i
    if(length(group)>1){
      warning("More than one group provided for ", i, ". First group being used.")
      group <- group[1]
    } 
    parameters[[i]][["group"]] <- group
  }

  group_names <- sapply(parameters[!names(parameters) %in% c("intercept","interactions")],function(x) x[["group"]])

  ## check data_structure
  if(is.null(data_structure)){
    if(any(!group_names %in% reserved_p_names)) stop("data_structure must be specified if there are more groups than 'intercept', 'observation', 'residual' and 'interactions' in parameter list", call.=FALSE)
  }else{
    if(!(is.matrix(data_structure)|is.data.frame(data_structure))) stop("data_structure is not a matrix or data.frame", call.=FALSE)   
  } 

  # User has to specify a "residual" level
  if(! "residual" %in% group_names) stop("One of the parameters groups must be 'residual'", call.=FALSE)
 
  # Check whether all groups match ones in data structure - If not, give error
  if(any(!group_names %in% c(colnames(data_structure),reserved_p_names))) stop("Group names in parameter list do not match group names in data_structure", call.=FALSE)

  # Check no group is called observation or residual - If not, give error
  if(any(colnames(data_structure) %in% reserved_p_names)) stop("'intercept', 'observation', 'residual' and 'interactions' are reserved names for grouping factors. Please rename grouping factors in data_structure", call.=FALSE)


########
##
########


  #i <- names(parameters)[1]



  for (i in names(parameters)[!names(parameters) %in% c("intercept","interactions")]){
    # p <- parameters[[i]]
    
    ## make cov from vcorr
    if(!is.null(parameters[[i]][["vcorr"]])){

      if(!is.null(parameters[[i]][["vcov"]])){
        message("vcov and vcorr are both specified for '",i,"', only vcov is being used")
      }else{
        if(!is.matrix(parameters[[i]][["vcorr"]])) stop("vcorr needs to be a matrix for ",i, call.=FALSE)
        if(nrow(parameters[[i]][["vcorr"]])!=ncol(parameters[[i]][["vcorr"]])) stop("need square vcorr matrix for ",i, call.=FALSE)
        if(!isSymmetric(parameters[[i]][["vcorr"]])) stop("vcorr matrix should be symmetric for ",i, call.=FALSE)
        if(!all(diag(parameters[[i]][["vcorr"]])>=0)){ stop("Variances for ",i," must all be >=0", call.=FALSE)}
        if(any(parameters[[i]][["vcorr"]][lower.tri(parameters[[i]][["vcorr"]])]< -1 | parameters[[i]][["vcorr"]][lower.tri(parameters[[i]][["vcorr"]])]>1)){ stop("Correlations for ",i," must all be between -1 and 1", call.=FALSE)}

        sd <- sqrt(diag(parameters[[i]][["vcorr"]]))
        corr <- parameters[[i]][["vcorr"]]
        diag(corr) <- 1
        parameters[[i]][["vcov"]] <- diag(sd) %*% corr %*% diag(sd)
      }
    }

    if(is.null(parameters[[i]][["beta"]]) & is.null(parameters[[i]][["vcov"]])){
      stop("'beta' or 'vcov' needs to be specified for ", i, call.=FALSE)
    }

    # If cov is not a matrix, make it one. Need to do this before working out k, as code below requires a matrix
    # if its a matrix check its square and symmetric
    # if its a vector, make its the diagonal of a square matrix
    # if neither give error
    if(!is.null(parameters[[i]][["vcov"]])){
      if(is.matrix(parameters[[i]][["vcov"]])){
        if(nrow(parameters[[i]][["vcov"]])!=ncol(parameters[[i]][["vcov"]])) stop("need square vcov matrix for ",i, call.=FALSE)
        if(!isSymmetric(parameters[[i]][["vcov"]])) stop("vcov matrix should be symmetric for ",i, call.=FALSE)
          #any(x[lower.tri(x)] != x[upper.tri(x)])

        ## remove rows and columns associated with 0 variance components to check positive definite
        # if(any(diag(parameters[[i]][["vcov"]])==0)){
        #   v <- parameters[[i]][["vcov"]]
        #   if(any(eigen(v[which(diag(v)!=0),which(diag(v)!=0)])$values<0))stop("vcov matrix should be positive definite for ",i, call.=FALSE)
        # }
        if(any(eigen(parameters[[i]][["vcov"]])$values<0))stop("vcov matrix should be positive semi-definite for ",i, call.=FALSE)
        if(!all(diag(parameters[[i]][["vcov"]])>=0)){ stop("Variances for ",i," must all be >=0", call.=FALSE)}
      }else if(is.vector(parameters[[i]][["vcov"]])){
        if(!all(parameters[[i]][["vcov"]]>=0)){ stop("Variances for ",i," must all be >=0", call.=FALSE)}
        parameters[[i]][["vcov"]] <- if(length(parameters[[i]][["vcov"]])==1) as.matrix(parameters[[i]][["vcov"]]) else diag(parameters[[i]][["vcov"]])
      }else{
        stop("vcov must be a symmetric square matrix or a vector", call.=FALSE)
      }
    }
  

## maybe dont have special rule for residual, because if in multivariate model want residuals at one level but not another, need to be able to specify beta=0 

    # if(i=="residual"){
    #   if(any(!names(parameters[[i]])) %in% c("mean", "vcov","vcorr")){
    #     message("Only mean and vcov/vcorr will be used for residual")
    #   }
    #   if(is.null(parameters[[i]][["vcov"]])){
    #   stop("'vcov' needs to be specified for ", i, call.=FALSE)
    # }
    #   #parameters[[i]][["beta"]] <- NULL


    #     param_names <- c("names", "group", "mean", "vcov", "vcorr", "beta", "n_level", "fixed", "n_response", "covariate")

    # }


    # If beta is not a matrix, make it one. good for working out k and for simulations, as code below requires a matrix
    if(!is.null(parameters[[i]][["beta"]])){
      if(is.vector(parameters[[i]][["beta"]])){
        parameters[[i]][["beta"]] <- matrix(parameters[[i]][["beta"]])
      }else if(!is.matrix(parameters[[i]][["beta"]])){stop("'beta' in ", i, " should be a vector or matrix", call.=FALSE)
      }
      if(n_response != ncol(parameters[[i]][["beta"]])){ 
        stop("Number of columns in beta is not the same as n_response for ",i, call.=FALSE)
      }
      ##
      ## remove any colnames on beta
      # colnames(parameters[[i]][["beta"]])<-NULL
    }else{ 
    ## if the number of responses and the size of cov are the same, and beta is not specified, then beta = I, as assuming that the user is simulating random effects
    ## if not the same then matrix of 1s
      if(n_response == ncol(parameters[[i]][["vcov"]])){ 
        parameters[[i]][["beta"]] <- diag(n_response)
      }else{
        parameters[[i]][["beta"]] <- matrix(1, nrow= ncol(parameters[[i]][["vcov"]]), ncol= n_response)  
      }
    }  
    ##
    # }else{ 
    #   if(is.null(parameters[[i]][["n_response"]])){ 
    #     parameters[[i]][["n_response"]] <- 1
    #   }else if(parameters[[i]][["n_response"]]>1){ parameters[[i]][["beta"]] <- diag(parameters[[i]][["n_response"]])
    #   }
    # } 

    # Work out number of variables at that level (k)
    # Check that size (k) of names, mean, cov, sd and var match - if not give error
    lengths <- c(length(parameters[[i]][["names"]]),
    	length(parameters[[i]][["mean"]]),
    	ncol(parameters[[i]][["vcov"]]),
    	nrow(parameters[[i]][["beta"]]), ## possibly change this if allowing matrix of sds for multivariate
      length(parameters[[i]][["functions"]])
    )
    k <- unique(lengths[lengths>0])
    if(length(k) != 1) stop("The number of parameters given for ", i, " are not consistent", call.=FALSE)
     

    # Check whether names specified
    # If not, generate names (length k)
    if(is.null(parameters[[i]][["names"]])){
      # if(k==1) parameters[[i]][["names"]] <- i
      # if(k>1) 
      parameters[[i]][["names"]] <- paste0(i,if(parameters[[i]][["group"]]!= "residual"){"_effect"},if(k>1){1:k})
    }else if(!is.vector(parameters[[i]][["names"]])){
      stop("'names' should be a vector for ", i, call.=FALSE)
    }
    # check everything is not just interaction terms
    # if(!any(!grepl(":",parameters[[i]][["names"]]))){
    #  stop("'names' only include interaction terms for ", i, call.=FALSE) 
    # }


    ## set n_level - assume that it is not input by user, if it is, it will be over-written
    parameters[[i]][["n_level"]] <- 
      if(parameters[[i]][["group"]] %in% c("observation","residual")){
        if(is.null(data_structure)){ 
          n
        }else{
          nrow(data_structure)
        }
      }else{
        length(unique(data_structure[,parameters[[i]][["group"]]]))
      }

    ## covariate and fixed should be fixed to false for observation and residual
    if(parameters[[i]][["group"]] %in% c("observation","residual")){
      parameters[[i]][["fixed"]] <- FALSE
      parameters[[i]][["covariate"]] <- FALSE
    }

    ## fixed    
    if(is.null(parameters[[i]][["fixed"]])) parameters[[i]][["fixed"]] <- FALSE
    
    if(parameters[[i]][["fixed"]] & parameters[[i]][["n_level"]] != k) stop("If fixed=TRUE, number of parameters should match the number of levels in grouping factor", call.=FALSE)
    
    if(parameters[[i]][["fixed"]] & is.null(parameters[[i]][["beta"]])) stop("If fixed =TRUE, beta also needs to be specified", call.=FALSE)

    if(parameters[[i]][["fixed"]] & (!is.null(parameters[[i]][["mean"]]) || !is.null(parameters[[i]][["vcov"]]) || !is.null(parameters[[i]][["functions"]]))) warning("Fixed=TRUE for ",i,", so mean, cov and functions are ignored", call.=FALSE)


    # Check whether covariate is specified
    if(is.null(parameters[[i]][["covariate"]])) parameters[[i]][["covariate"]] <- FALSE
  
    if(parameters[[i]][["covariate"]] & (!is.null(parameters[[i]][["mean"]]) || !is.null(parameters[[i]][["vcov"]]) || !is.null(parameters[[i]][["functions"]]))) warning("Covariate=TRUE for ",i,", so mean, cov and functions are ignored", call.=FALSE)
    
    if(parameters[[i]][["covariate"]] && is.null(parameters[[i]][["beta"]])) stop("If covariate =TRUE, beta also needs to be specified", call.=FALSE)

    if(parameters[[i]][["covariate"]] && is.null(parameters[[i]][["fixed"]])) stop("covariate =TRUE and fixed=TRUE for ", i, call.=FALSE)

    if(parameters[[i]][["covariate"]] && !is.numeric(data_structure[,parameters[[i]][["group"]]])) stop("If covariate =TRUE, the corresponding grouping factor in the data_structure needs to be coded as a numeric", call.=FALSE)

    # Check whether mean specified
    # If not, rep(0,k)
    if(is.null(parameters[[i]][["mean"]])){
      parameters[[i]][["mean"]] <- rep(0,k)
    }else if(!is.vector(parameters[[i]][["mean"]])){
      stop("'mean' should be a vector", call.=FALSE)
    }
    
    # Check whether cov specified
    # If not, diag(k)
    if(is.null(parameters[[i]][["vcov"]])) parameters[[i]][["vcov"]] <- diag(k)

    
    ### functions
    if(is.null(parameters[[i]][["functions"]])){
      parameters[[i]][["functions"]] <- rep("identity",k)
    }else if(!is.vector(parameters[[i]][["functions"]])){
      stop("'functions' should be a vector for ", i, call.=FALSE)
    }else if(!all(is.character(parameters[[i]][["functions"]]))){
      stop("'functions' should be a character vector for ", i, call.=FALSE)
    }else{
      parameters[[i]][["functions"]][is.na(parameters[[i]][["functions"]])] <- "identity"
    }
    ## should they be input as characters?

  ## add names to means, cov and betas
    names(parameters[[i]][["functions"]]) <- names(parameters[[i]][["mean"]]) <- rownames(parameters[[i]][["beta"]]) <- rownames(parameters[[i]][["vcov"]]) <- colnames(parameters[[i]][["vcov"]]) <- parameters[[i]][["names"]]

  }

  ##check whether all betas have same dimension
  # j <- n_phenotypes(parameters)


########
## intercept
########
  if(!is.null(parameters[["intercept"]])){
    ## needs to be same length as n_response
    if(length(parameters[["intercept"]])!=n_response) stop("intercept needs to same length as n_response")
    if(!is.vector(parameters[["intercept"]])) stop("intercept needs to be a vector of length N_reponse") 

  }else{
    parameters[["intercept"]] <- rep(0, n_response)
  }


########
## interactions
########

  pred_names <- c(do.call(c, lapply(parameters[!names(parameters) %in% c("intercept")],function(x) x[["names"]])), colnames(known_predictors[["predictors"]]))
  if(any(duplicated(pred_names))) stop("Predictor names must be unique", call.=FALSE)
  
# pred_names <- c("a","b","c","d")
# interactions <- list(names=c("a:b","c:d"), beta=c(2,3), cov=2)

  if(!is.null(parameters[["interactions"]])){
    interactions <- parameters[["interactions"]]

    if(is.null(interactions[["names"]])) stop("'names' need to be specified for interactions", call.=FALSE) 
    if(any(!names(interactions) %in% c("names", "beta"))){
        message("Only 'names' and 'beta' will be used for interactions")
        interactions <- interactions[c("names", "beta")]
      }
    interaction_names <- interactions[["names"]]
    interaction_variables <- unique(c(strsplit(interaction_names,":"), recursive=TRUE))

    ## check that all main effects are also specified along with interactions
    if(!all(interaction_variables %in% pred_names )) stop("variables included in interactions are not all specified in the 'names' arguments of the parameter list", call.=FALSE) 


    if(!is.null(interactions[["beta"]])){
      if(is.vector(interactions[["beta"]])){
        interactions[["beta"]] <- matrix(interactions[["beta"]])
      }else if(!is.matrix(interactions[["beta"]])){stop("'beta' in interactions should be a vector or matrix", call.=FALSE)
      }
      if(n_response != ncol(interactions[["beta"]])){ 
        stop("number of columns in beta is not the same as n_response for interactions", call.=FALSE)
      }
    }else{ 
      interactions[["beta"]] <- matrix(1,length(interaction_names),n_response)
      #diag(interactions[["n_response"]])
    } 

    ## names and betas needs to be vectors of same dimensions
    if(length(interactions[["names"]])!=nrow(interactions[["beta"]])) stop("'names' and 'beta' need to be the same length in interactions ", call.=FALSE) 
    
    rownames(interactions[["beta"]]) <- interaction_names
    parameters[["interactions"]][["beta"]] <- interactions[["beta"]]
  }







  ##Check extra parameters
  param_names <- c("names", "group", "mean", "vcov", "vcorr", "beta", "n_level", "fixed", "covariate", "functions")

  e_p <- unlist(sapply(parameters, function(x){
    names(x)[!names(x) %in% param_names]
    }))

  if(length(e_p)>0){
    ## check is extra parameters are vectors
    e_p_vector <- !unlist(sapply(parameters, function(x){
      sapply(x[!names(x) %in% param_names],is.vector)
      }))
    if(any(e_p_vector)) stop("Additional parameters given to parameters lists must be vectors, this is not the case for ",names(e_p_vector)[e_p_vector], call.=FALSE)
  
    ## check length of all extra parameters is 1
    e_p_length <- unlist(sapply(parameters, function(x){
      sapply(x[!names(x) %in% param_names],length)
      }))
    if(any(e_p_length>1)) stop("Additional parameters given to parameters lists must be length 1, this is not the case for ",names(e_p_length)[e_p_length>1], call.=FALSE)

  }


  ## Check whether all names in data_structure and parameters contain only words, digits, : and _
  all_names <- c(response_names, e_p, pred_names, colnames(data_structure))

  if(any(duplicated(all_names))) stop("Response names, names in data structure and names in parameters must be unique")

  if(!all(grepl("^[A-z0-9_:]*$",all_names))) stop("Response names, names in data structure and names in parameters must be alphanumeric, '_' or ':'", call.=FALSE)

  ### check no names are repeated!!
  if(any(duplicated(e_p)) || any(e_p %in% c(pred_names, colnames(data_structure)))) stop("Additional parameters names must be unique", call.=FALSE)

  ### check all names have at least 1 character!!
  if(any(nchar(c(pred_names, e_p, colnames(data_structure)))==0 )) stop("Specified names must have nchar>0", call.=FALSE)


	return(parameters)

}






# data(BTdata)
# known_predictors <- list(predictors=BTdata[,c("hatchdate","tarsus")], beta=matrix(c(1,2,3),nrow=3))
# n_response <- 1

### Function to checks known_predictors list

fill_preds <- function(known_predictors,n, n_response,...){
  
  if(!is.list(known_predictors)) stop("known_predictors must be a list", call.=FALSE)


  ## checks on predictors
  preds <- known_predictors[["predictors"]]
  if(is.null(preds)) stop("known_predictors must contain predictors ", call.=FALSE)  
  if(!(is.data.frame(preds) | is.matrix(preds)))stop("predictors in known_predictors must be a data.frame or matrix", call.=FALSE)
  if(n!=nrow(preds)) stop("The number of observation specified does not match the number of rows in known_predictors", call.=FALSE)
  if(!length(colnames(preds))>0) stop("predictors in known_predictors must have column names", call.=FALSE)
  if(is.data.frame(preds)) preds <- as.matrix(preds)
  if(any(!apply(preds,2,is.numeric))) stop("predictors in known_predictors must be numeric", call.=FALSE)
  
  ## checks on predictors
  betas <- known_predictors[["beta"]]
  if(!is.null(betas)){
    if(is.vector(betas)){
      betas <- matrix(betas)
    }else if(!is.matrix(betas)) {
      stop("'beta' in known_predictors should be a vector or matrix", call.=FALSE)
    }
    if(!is.numeric(betas))stop("beta in known_predictors must be numeric", call.=FALSE)
    if(n_response != ncol(betas)) stop("Number of columns in beta is not the same as n_response for known_predictors", call.=FALSE)      
    if(ncol(preds) != nrow(betas)) stop("The dimensions of beta do not match the number of predictors for known_predictors", call.=FALSE)
      ##
  }else{ 
    betas <- matrix(1, nrow= ncol(preds), ncol= n_response)  
  }

  rownames(betas) <- colnames(preds)

  return(list(predictors=preds,beta=betas))
}
# kn<-fill_preds(list(predictors=BTdata[,c("hatchdate","tarsus")], beta=matrix(c(1,2,3),ncol=1)),1)

