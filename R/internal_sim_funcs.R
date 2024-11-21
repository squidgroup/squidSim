### wrapper for MCMCglmm::rbv that allows 0 variances 
rbv0 <- function(pedigree, G){
  X <- matrix(0, nrow=nrow(pedigree), ncol=nrow(G))
  index <- which(diag(G)!=0)
  if(length(index)>0){
    if(any(diag(G)==0)) G <- G[index,index]
    X2 <- MCMCglmm::rbv(pedigree=pedigree, G=G)
    X[,index] <- X2 
  }
  X
}

### wrapper for mvnfast::rmvn that allows 0 variances 
rmvn0 <- function(n,mu,sigma){
  X <- matrix(0, nrow=n, ncol=nrow(sigma))
  index <- which(diag(sigma)!=0)
  X[,index] <- mvnfast::rmvn(n=n, mu=mu[index], sigma=sigma[index,index])
  X
}


### function to generate ar1 matrix 
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
      (1:n - 1))
  rho^exponent
}
# from https://www.r-bloggers.com/2020/02/generating-correlation-matrix-for-ar1-model/
#ar1_cor(n=10,rho=0.5)


## cov_str_check
##  performs error checking of all cov structures
cov_str_check <- function(data_structure, pedigree, phylogeny, cov_str, parameters,...){
  cs_check <- lapply(c("pedigree","phylogeny","cov_str"), function(j){
      cs <- get(j)
      if(!is.list(cs) | is.data.frame(cs)) stop(j, " needs to be a list", call.=FALSE) 
    })

  param_names <- names(parameters)
  group_names <- sapply(parameters,function(x) x$group)
  ped_names <- names(pedigree)
  phylo_names <- names(phylogeny)
  cov_names <- names(cov_str)
  cs_names <- c(ped_names,phylo_names,cov_names)

  if(any(duplicated(cs_names))) stop("Cannot have multiple covariance structures (pedigree/phylogeny/cov_str) linking to the same item in the parameter list. If multiple covariance structures are needed to be linked to the same grouping factor in the data_structure (for example simulating additive genetic and dominance effects), then create multiple items in the parameter list, with different names, but the same 'group', and link the covariance structures accordingly.", call.=FALSE)
  if(any(!cs_names %in% param_names)) stop("Some names in pedigree/phylogeny/cov_str are not in the parameter list", call.=FALSE)


  lapply(colnames(data_structure), function(i){
      # i = colnames(data_structure)[1]
    
    # are any of the parameter list names associated with it     
    list_names <- param_names[group_names %in% i]
    
    ped_link <- list_names[list_names %in% ped_names]
    phylo_link <- list_names[list_names %in% phylo_names]
    cov_link <- list_names[list_names %in% cov_names]
    all_link <- c(ped_link,phylo_link,cov_link)

    if(length(all_link)>1){ 
      warning("Multiple covariance structures linked to ",i,". The function assumes that the covariance structures are ordered exactly the same. If they are not then the simulations will not run as you expect. You will need to create multiple columns in grouping structure and link different covariance structures to each one.", call.=FALSE) 
    }
  })

  # ped_names <- c(names(pedigree),names(phylogeny),names(cov_str))

  # names_check <- lapply(ped_names,function(i){
  # # data_structure[,i]
  #   if(!all(unique(rownames(chol_str_all[[i]])) %in% unique(data_structure[,parameters[[i]]$group]))) stop(paste("all IDs in the pedigree/phylogeny/cov_str linked with", i, "are not in the data_structure"), call.=FALSE)
  #   if(!all(unique(data_structure[,parameters[[i]]$group]) %in% unique(rownames(chol_str_all[[i]])))) stop(paste("all IDs in data_structure are not in the pedigree/phylogeny/cov_str linked with", i), call.=FALSE)      
  # })
  
}



### function to turn data_structure into indexes. Matches names in data_structure to linked pedigree/phylogeny/cov_str to make sure indexing is correct. Doesnt do any error checking
index_factors <- function(data_structure, pedigree, phylogeny, cov_str, parameters, index_link,suppress_index_warning,...){
  
  p_names <- names(parameters)[!names(parameters)%in%c("intercept","interactions")]

  group_names <- sapply(p_names,function(x) parameters[[x]]$group)

  new_ds <- if(is.null(data_structure)){   
    
    NULL

  }else if(is.null(pedigree) & is.null(phylogeny) & is.null(cov_str)){

    apply(data_structure,2,function(x) as.numeric(factor(x)))

  }else{

    do.call(cbind,
    ## for each column of the data structure
      lapply(colnames(data_structure), function(i){
      # i = colnames(data_structure)[1]
        

        # are any of the parameter list names associated with it     
        list_names <- group_names[group_names %in% i]
        
        ped_link <- list_names[list_names %in% names(pedigree)]
        phylo_link <- list_names[list_names %in% names(phylogeny)]
        cov_link <- list_names[list_names %in% names(cov_str)]
       
        # if linked with a something get the row number in the relevant cov_str, so the indexing will match the order in that cov str.
        # Can only match with one so is assuming that if multiple things are linked with it that the ordering is the same, so takes the first linked cov str it can find and indexes according to that

        if(length(ped_link)>0){
          # first column of pedigree
          match(data_structure[,i],pedigree[[ped_link[1]]][,1])
        }else if(length(phylo_link)>0){ 
          # names of phylogeny
          match(data_structure[,i],phylogeny[[phylo_link[1]]]$tip.label)
        }else if(length(cov_link)>0){
          # rownames(of cov_str) 
          match(data_structure[,i],rownames(cov_str[[cov_link[1]]]))
        }else{
          as.numeric(factor(data_structure[,i]))
        }
      })
    )
    
  }
  colnames(new_ds) <- colnames(data_structure)

  if(!is.null(index_link)){
    more_ds <- do.call(cbind,lapply(index_link, function(x){
      ds_links <- strsplit(x, '-')[[1]]
      
      ## give warning if not all levels match

      if(!all(data_structure[,ds_links[1]] %in% data_structure[,ds_links[2]]) & !suppress_index_warning) warning(paste("Not all levels are of", ds_links[1], "are present in", ds_links[2], "meaning that there will be NAs in the new grouping factor"), call.=FALSE)
      
      new_ds[,ds_links[2]][match(data_structure[,ds_links[1]], data_structure[,ds_links[2]])]

      }
    ))
    return(cbind(new_ds,more_ds))
  }else{
    return(new_ds)
  }

}


cov_str_list <- function(parameters, data_structure, phylogeny, phylogeny_type, cov_str,...){
#pedigree, pedigree_type, 

  # ped_chol <- sapply(names(pedigree), function(x){
  #   if(pedigree_type[[x]]=="A") Matrix::chol(nadiv::makeA(pedigree[[x]]))
  #   else if(pedigree_type[[x]]=="D") Matrix::chol(nadiv::makeD(pedigree[[x]]))
  #   else if(pedigree_type[[x]]=="E") Matrix::chol(nadiv::makeAA(pedigree[[x]]))
  # })

  phylo_chol <- sapply(names(phylogeny), function(x){
    phylo_vcv <- ape::vcv(phylogeny[[x]], corr = TRUE)
    # if(phylogeny_type[[x]]=="OU") {
    #   ## need way of specifying alpha
    #   phylo_vcv <- exp(phylo_vcv * - alpha)
    #   diag(phylo_vcv) <- 1
    # }
    methods::as(chol(phylo_vcv), "dgCMatrix")
  })

  cor_chol <- lapply(cov_str, function(x) methods::as(chol(x), "dgCMatrix"))

  chol_str_all<-c(phylo_chol,cor_chol)#ped_chol,

  add_list<-names(parameters)[!names(parameters) %in% c(names(chol_str_all),"intercept","interactions")]
  for(i in add_list){
    chol_str_all[[i]] <- NULL#Matrix::Diagonal(parameters[[i]][["n_level"]])
  }
  return( chol_str_all)
}



sim_predictors <- function(parameters, str_index, cov_str_all, known_predictors, pedigree, ...){

  traits <- do.call(cbind, lapply( names(parameters)[!names(parameters)%in%c("intercept","interactions")], function(i){  

# i<-"animal"
    p <- parameters[[i]]
    k <- length(p$mean)
    n <- p$n_level

    if(p$fixed){
      
      ## if factor make design matrix
      # x<-stats::model.matrix(stats::formula(paste("~ factor(",p$group,")-1")),as.data.frame(str_index))
      ## if there are group names in the data structure, then use factor levels with the same order as the specified names
      ds_levels_group <- unique(data_structure[,p[["group"]]])

      if(all(ds_levels_group %in% 1:length(ds_levels_group)) ){
        fac_levels <- 1:length(ds_levels_group)
      }else{
        fac_levels <- p$names
      }

      x<-stats::model.matrix(stats::formula(paste("~ factor(",p$group,",levels=fac_levels)-1")),as.data.frame(data_structure))

    }else if(p$covariate){
      
      ## if covariate, make design matrix from data structure
      x<- matrix(rep(str_index[,p$group],k),nrow(str_index),k)

    }else{
      ## otherwise simulate 'traits' at each level from multivariate normal 

      # x <- if(is.null(cov_str_all[[i]])) {
      #   mvnfast::rmvn(n=n, mu=p$mean, sigma=p$vcov)
      # }else{
      #   methods::as(Matrix::crossprod(cov_str_all[[i]],mvnfast::rmvn(n=n, mu=p$mean, sigma=p$vcov)),"matrix")
      # }

      if(i %in% names(pedigree)){
        ## if name is listed in pedigree argument, link to pedigree
        x <- rbv0(pedigree[[i]],p$vcov)
      }else{
        x <- rmvn0(n=n, mu=p$mean, sigma=p$vcov)  
      }
      
      if(!is.null(cov_str_all[[i]])) {
        x <- methods::as(Matrix::crossprod(cov_str_all[[i]],x),"matrix")
      }

      ### apply functions
      ### add them in likes betas, making sure they are in the right order? then apply them to the cols?
      x <- sapply(1:k,function(i) sapply(x[,i],p[["functions"]][i]))
      
      ## expand traits to be the same length as the number of observations using data structure  
      if(!p$group %in% c("observation","residual")) x <- x[str_index[,p$group],,drop=FALSE]
    }

    ## use names form parameter list 
    colnames(x) <- p$names

    return(x)
  }))


# traits <- data.frame(a=rnorm(100),b=rnorm(100),c=rnorm(100),d=rnorm(100))
# interactions <- list(names=c("a:b","c:d", "exp(a)"), beta=c(2,3))

  ## add in known predictors
  if(!is.null(known_predictors)){
    traits<-cbind(traits,known_predictors$predictors)
  }

  ## add in interactions
  if(!is.null(parameters[["interactions"]])){

    x_int <- do.call(cbind,lapply(strsplit(parameters[["interactions"]]$names,":"), function(j){
        eval(parse(text=paste(j, collapse="*")), envir = as.data.frame(traits) )
    }))

    colnames(x_int) <- parameters[["interactions"]]$names

    traits<-cbind(traits,x_int)

  }
  

  return(traits)

}




generate_y <- function(predictors, intercepts, betas, str_index, model, y_pred_names,extra_param,...){

  ## evaluate model
  ## - if model is missing, add all simulated predictors together
  if(is.null(model)) {
    y <- predictors %*% betas
    y <- t(t(y) + intercepts)
  } else {
    ## for evaluation with model formula 
    
    beta_predictors <- predictors * rep(betas, rep(nrow(predictors),length(betas))) #predictors %*% diag(as.vector(betas))
    y_predictors <- cbind(beta_predictors,predictors,str_index)
    colnames(y_predictors) <- y_pred_names

    ## allow I() and subsets to be properly linked to y_predictors
    model <- gsub("I\\((\\w+)\\)","\\1_raw",model)
    model <- gsub("\\[(\\w+)\\]","\\[\\1_ID\\]",model)

    # evaluate the formula in the context of y_predictors and the extra params
    # y <- eval(parse(text=model), envir = c(as.data.frame(y_predictors),as.list(extra_param)))
    
    model2 <- paste(model,";\n return(data.frame(mget(ls()[!ls() %in% c(colnames(y_predictors),names(extra_param),'intercept')])))")
    y <- eval(parse(text=model2), envir = c(as.data.frame(y_predictors),intercept=intercepts,as.list(extra_param)))

    if(!grepl("intercept",model)) y <- y + intercepts
    
    if(is.vector(y)) y <- matrix(y)
  }
  
  return(y)
}

generate_y_list <- function(parameters, str_index, predictors, model,known_predictors,...){
  
  if(missing(str_index)) str_index <- NULL

  intercepts <- parameters[["intercept"]]
  parameters <- parameters[!names(parameters) %in% c("intercept")]
  
  ## put all betas together
  #order betas to match predictors
  betas <- do.call(rbind,lapply(parameters,function(x) x$beta))
  if(!is.null(known_predictors)){
    betas<-rbind(betas,known_predictors$beta)
  }
  betas<-betas[colnames(predictors[[1]]),]

  y_pred_names <- c(colnames(predictors[[1]]), paste0(colnames(predictors[[1]]),"_raw"), if(!is.null(str_index)){paste0(colnames(str_index),"_ID")})

  ## extract and name extra parameters
  if(!is.null(model)){

    param_names <- c("names", "group", "mean", "vcov", "vcorr", "beta", "fixed", "covariate", "n_level", "functions")
    extra_param <- unlist(sapply(parameters, function(x){ x[!names(x) %in% param_names] }))
        
    if(!is.null(extra_param)){
      names(extra_param) <- unlist(sapply(parameters, function(x) names(x)[!names(x) %in% param_names]
        ))
      ## check extra param names dont clash with y_trait names
      if(any(names(extra_param) %in% colnames(y_pred_names))) stop("You cannot name extra parameters the same as any variables")
    }
  }

  y <- lapply(predictors, function(x) generate_y(x, intercepts=intercepts, betas=betas, str_index=str_index,  model=model, y_pred_names=y_pred_names,extra_param=extra_param))

  return(y)

}


transform_dist <- function(y, family, link, response_names,...){

  inv <- function(x) 1/x

  j <- ncol(y)

  if(length(link)==1 & j>1) link <- rep(link,j)
  if(length(family)==1 & j>1) family <- rep(family,j)

  ## convert the link argument into an actual function
  link_function <- 
  ifelse(link=="log", "exp", 
  ifelse(link=="inverse", "inv", 
  ifelse (link=="logit", "plogis", 
  ifelse (link=="probit", "pnorm", 
   link))))
  

  y_family <-  sapply(1:j,function(i){
    ## apply link function to y
  y_link <- get(link_function[i])(y[,i])
    ## sample from poisson or binomial 
    if(family[i]=="gaussian") y_link else 
    if(family[i]=="poisson") stats::rpois(length(y_link),y_link) else 
    if(family[i]=="binomial") stats::rbinom(length(y_link),1,y_link)
  })
  
  ## think colnames(y) will be null unless model is specified but should check :)
  if(is.null(colnames(y))){
    if(!is.null(response_names)){
      colnames(y_family) <- response_names
    }else{
      colnames(y_family) <- if(j==1)"y" else paste0("y",1:j)  
    }
  }else{
    colnames(y_family) <- colnames(y)
  }

  return(y_family)
}

