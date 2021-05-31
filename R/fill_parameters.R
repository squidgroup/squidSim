
## I've used loops rather than apply functions in here because then the original parameter list can then be added to rather than new lists made - this will be slightly slower but very negligible given their size
fill_parameters <- function(parameters,data_structure){

  # Check whether list is given
  if(!is.list(parameters)) stop("parameters are not provided as a list")

  	# if(exists("data_structure"))

  ## check data_structure
  if(!(is.matrix(data_structure)|is.data.frame(data_structure))) stop("data_structure is not a matrix or data.frame")

  # Check whether group specified -If not, then give name in list
  # Check whether length(group)==1 - If not, give warning and use only first group

  #i <- names(parameters)[1]
  for (i in names(parameters)){
    group <-parameters[[i]]$group
    if(is.null(group)) group <- i
    if(length(group)>1){
      warning("More than one group provided for ", i, ". First group being used.")
      group <- group[1]
    } 
    parameters[[i]]$group <- group
  }

  # User has to specify a "residual" level
  if(! "residual" %in% sapply(parameters,function(x) x$group)) stop("One of the parameters groups must be 'residual'")
 
  # Check whether all groups match ones in data structure - If not, give error
  if(any(!sapply(parameters,function(x) x$group) %in% c(colnames(data_structure),"residual"))) stop("Group names in parameter list do not match group names in data_structure")

#i <- names(parameters)[1]
  for (i in names(parameters)){
    
  # Work out number of variables at that level (k)
  # Check that size (k) of names, mean, cov, sd and var match - if not give error

    lengths <- c(length(parameters[[i]]$names),
    	length(parameters[[i]]$mean),
    	ncol(parameters[[i]]$cov),
    	length(parameters[[i]]$sd), ## possibly change this if allowing matrix of sds for multivariate
    	length(parameters[[i]]$var)
    )
    k <- unique(lengths[lengths>0])
    if(length(k) != 1) stop("The number of parameters given for ", i, " are not consistent")
    
    
    # Check whether names specified
    # If not, generate names (length k)
    if(is.null(parameters[[i]]$names)){
      if(k==1) parameters[[i]]$names <- i
      if(k>1) parameters[[i]]$names <- paste(i,1:k,sep="_")
    }

    # Check whether mean specified
    # If not, rep(0,k)
    if(is.null(parameters[[i]]$mean)) parameters[[i]]$mean <- rep(0,k)
    
    # Check whether cov specified
    # If not, diag(k)
    if(is.null(parameters[[i]]$cov)) parameters[[i]]$cov <- diag(k)

    # Check whether sd and var specified - If both, give error
    if(!is.null(parameters[[i]]$sd) & !is.null(parameters[[i]]$var)) stop("Specify either sd or var, not both")    
    # - If neither, sd=rep(1,k)
    if(is.null(parameters[[i]]$sd) & is.null(parameters[[i]]$var)) parameters[[i]]$sd <- rep(1,k)
    # - If var, sd=sqrt(var)
    if(is.null(parameters[[i]]$sd)) parameters[[i]]$sd <- sqrt(parameters[[i]]$var)

    ## Check whether number of levels is specified
    # - if no take from data structure 
    # - if yes check it matches data structure 	
    if(is.null(parameters[[i]]$n_level)){
	  if(i=="residual") {
          parameters[[i]]$n_level <- nrow(data_structure)
      } else {
		  parameters[[i]]$n_level <- length(unique(data_structure[,parameters[[i]]$group]))
      }
	} else {
      
	}

  
  }

	return(parameters)



}



### possibly make it so that if data_structure is not specified then make one with completely crossed random effects?
### then you would need to specify the number of levels