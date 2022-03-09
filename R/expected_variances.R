## function to compute expected variance of the product of variables. if two variables are given, accounts for non-independence of variables (i.e. non-0 covariance), but for more than 3 assumes that they are independent (don't' think there is a generalisation for the produce of 3 or more dependent variables)
prod_var <- function(means,vcov){
	vars <- diag(vcov)
	
	if(length(means)>2){
		prod(vars + means^2) - prod(means)^2
		# https://stats.stackexchange.com/questions/52646/variance-of-product-of-multiple-independent-random-variables
	}else{
		c <- vcov[lower.tri(vcov)]
		prod(vars+means^2) + c^2 + 2*c*prod(means) - prod(means)^2
		# https://stats.stackexchange.com/questions/15978/variance-of-product-of-dependent-variables
		# https://math.stackexchange.com/questions/1889402/covariance-of-two-squared-not-zero-mean-random-variables
	}
	## this is an equation for variance of k dependent random variables centered on 0 - probably two complex to worry about
	## https://stats.stackexchange.com/questions/60414/variance-of-product-of-k-correlated-random-variables
}

prod_means <- function(means,vcov){
	##https://www.physicsforums.com/threads/expectations-on-the-product-of-two-dependent-random-variables.276125/
	
	if(length(means)>2){
		prod(means)
	}else{
		prod(means) + vcov[lower.tri(vcov)]
	}
}

make_big_matrix<-function(x){
	all_names <- c(sapply(x, function(i) colnames(i) ), recursive=TRUE)
	mat_index <- c(0,cumsum(sapply(x, function(i) nrow(i))))
	
	mat <- diag(max(mat_index))
	colnames(mat) <- rownames(mat) <-all_names

	for(i in 2:length(mat_index)){
		indexes<-(mat_index[i-1]+1):mat_index[i]
		mat[indexes,indexes] <- x[[i-1]]
	}
	mat
}

#' @title expected_variance
#' @description Calculate expected mean and variance in response variable(s)
#'
#' @param squid A squid object created with the simulate_population() function
#' @details 
#' Calculates the expect variance from the simulation parameters,. as well as breaking down the source of this variance into the hierarchical levels given in the parameter list, and into individual predictors. This function has several limitations. 
#' 1. Doesn't work when 'model' is specified in simulate_population()
#' 2. Assumes grouping factors in the data_structure are balanced - unbalanced designs can change the observed variance
#' 3. Will be inaccurate with transformed variables (i.e. if functions are specified)
#' 4. Doesn't deal well with three way interactions and over, unless they are the same variable (i.e. polynomials)
#' 5. Doesn't account for covariance between interaction terms, e.g. if there is are two interaction terms rain:temp and wind:temp, they will covary, but this additional variance in the response is not calculated
#' @author Joel Pick - joel.l.pick@gmail.com
#' @return A list with means and variance at each level
#' @examples
#' \dontrun{
#' squid_data <- simulate_population(
#' data_structure = make_structure(structure = "individual(10)", repeat_obs=2),
#'	parameters=list(
#'    individual=list(
#'  		vcov=1.2
#'  	),
#'  	observation=list(
#'      names=c("temperature","rainfall", "wind"),
#'      mean = c(10,1 ,20),
#'      vcov =matrix(c(
#'        1, 0, 1,
#'        0,0.1,0,
#'        1, 0, 2
#'        ), nrow=3 ,ncol=3),
#'      beta =c(0.5,-3,0.4)
#'    ),
#'    residual=list(
#'      mean=10,
#'      vcov=1
#'    )
#'  )
#')
#' 
#' expected_variance(squid_data)
#' }
#' @export


expected_variance <- function(squid){
	intercept <- squid$param$intercept
	param <- squid$param[names(squid$param)!="intercept"]
	data_structure <- squid$data_structure
	known_predictors <- squid$known_predictors

	if(any(sapply(param, function(i) any(i$functions!="identity")))){
		message("This will be inaccurate with transformed variables (i.e. using the functions argument)")
	}

	p_names <- names(param)[names(param)!="interactions"]
	
	fixed <- p_names[sapply(p_names, function(i) param[[i]]$fixed)]
	if(length(fixed)>0){
		for (i in fixed){
			levels<-as.vector(table(data_structure[,param[[i]]$group]))
			p <- levels/sum(levels)
			param[[i]]$mean <- p
			param[[i]]$vcov <- diag(p*(1-p), nrow=length(levels))
			names(param[[i]]$mean) <- colnames(param[[i]]$vcov) <- rownames(param[[i]]$vcov) <- param[[i]]$names
		}
	}
	covariate <- p_names[sapply(p_names, function(i) param[[i]]$covariate)]
	if(length(covariate)>0){
		for (i in covariate){
			z <- data_structure[,param[[i]]$group]
			param[[i]]$mean <- rep(mean(z), length(param[[i]]$names))
			param[[i]]$vcov <- diag(stats::var(z), length(param[[i]]$names))
			names(param[[i]]$mean) <- colnames(param[[i]]$vcov) <- rownames(param[[i]]$vcov) <- param[[i]]$names
		}
	}

	if(!is.null(known_predictors)){
		param$known_predictors <- list(
			mean = colMeans(known_predictors[["predictors"]]),
			vcov = cov(known_predictors[["predictors"]]),
			beta = known_predictors[["beta"]]
		)
	}

	if("interactions" %in% names(param)){
		
		means1 <- do.call(c,c(sapply(p_names, function(i) param[[i]]$mean ), use.names=FALSE))

		covs1 <- make_big_matrix(lapply(p_names, function(i) param[[i]]$vcov ))

		## if interaction names are the same, then cov = var
		## https://stats.stackexchange.com/questions/53380/variance-of-powers-of-a-random-variable
		## Var(𝑋𝑛)=𝔼[𝑋2𝑛]−𝔼[𝑋𝑛]2 
		int_names <- param[["interactions"]]$names
		int_var <- lapply(strsplit(int_names,":"), function(j){
			#For two way interactions, the expected means and variances take into account 
			if(length(j)>2) warning("For three way interactions and above, covariance between variables is ignored when calculating expected means and variances (with the exception of polynomials)")
			list(
				cov = prod_var(means1[j],covs1[j,j]),
				mean = prod_means(means1[j],covs1[j,j])
				)
		})

		param[["interactions"]]$vcov <- diag(sapply(int_var,function(x)x$cov), nrow=length(int_names))
		param[["interactions"]]$mean <- sapply(int_var,function(x)x$mean)
		names(param[["interactions"]]$mean) <- colnames(param[["interactions"]]$vcov) <- rownames(param[["interactions"]]$vcov) <- int_names
	}



	means <- do.call(c,c(lapply(param, function(p) p$mean ), use.names=FALSE))
	covs <- make_big_matrix(lapply(param, function(p) p$vcov ))
	betas <- do.call(rbind,lapply(param, function(p) p$beta ))

	total_var <- as.vector(t(betas) %*% covs %*% betas)
		
	out <- list( 
		## total
		total = c(
			mean = intercept + sum(betas * means),
			var = total_var
		),
		
		##hierarchy
		groups = data.frame(
			mean= c(intercept=intercept,sapply(param, function(p) sum(p$beta*p$mean))),
			var=c(intercept=0,sapply(param, function(p) t(p$beta) %*% p$vcov %*% p$beta ))
			),
	
		variables = data.frame(
			mean = rbind(intercept=intercept,betas * means),
			var = rbind(intercept=0,betas * covs %*% betas)
		)
	 	
	)
	
	class(out) <- "squid_var"
	return(out)
}


#' @title print.squid_var
#' @description Print method for class 'squid_var'
#' @param x an R object of class 'squid_var'
#' @param ... further arguments passed to or from other methods.
#' @export
print.squid_var <- function(x, ...){
	cat("Contribution of the simulated predictors to the mean and variance in the response\n\n")

	cat("Expected Mean:",x$total["mean"],"\n")
	cat("Expected Variance:",x$total["var"],"\n\n")

	cat("Contribution of different hierarchical levels to grand mean and variance:\n")
	print(x$groups)
	cat("\n\nContribution of different predictors to grand mean and variance:\n")
	print(x$variables)
}

