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
#' @description Calculate expected variance in response variable(s)
#'
#' @param squid A squid object
#' @details 
#' Calculated the expect variance from the simulation parameters. Has several limitaitons. 
#' 1. Doesn't work when 'model' is specific in simulate_population
#' 2. Assumes random factors are balanced - unbalanced designs will change the observed variance
#' 3. Will be inaccurate with transformed variables (i.e. if functions are specified)
#' 4. Doesn't deal well with three way interactions and over, unless they are the same variable (i.e. polynomials)
#' 5. Doesn't account for covariance between interaction terms, e.g. if there is are two interaction terms rain:temp and wind:temp, they will covary, but this additional variance in the response is not calculated
#' @author Joel Pick - joel.l.pick@gmail.com
#' @return A data.frame with the data structure
#' @examples
#' # simple data structure with 5 'individuals' and 2 observations per individual
#' make_structure(structure="individual(5)", repeat_obs=2)
#' 
#' # nested data structure with 2 sexes, 5 individuals per sex
#' # and 2 observations per individual
#' make_structure(structure="sex(2)/individual(5)", repeat_obs=2)
#' 
#' # crossed data structure with 5 individuals in 2 treatments 
#' # and 2 observations per individual and treatment combination
#' make_structure(structure="treatment(2) + individual(5)", repeat_obs=1)
#' 
#' @export


expected_variance <- function(squid){
	param <- squid$param
	p_names <- names(param)[names(param)!="interactions"]

	if(any(sapply(param, function(i) any(i$functions!="identity")))){
		message("This will be inaccurate with transformed variables (i.e. using the functions argument)")
	}
	
	if("interactions" %in% names(param)){
		
		p_names <- names(param)[names(param)!="interactions"]

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
list( 
	## total
	total = cbind(
		mean = sum(betas * means),
		var = total_var
		),
	
	##hierarchy
	groups = cbind(
		mean= sapply(param, function(p) sum(p$mean)),
		var=sapply(param, function(p) t(p$beta) %*% p$vcov %*% p$beta )),

	variables = cbind(
		mean = betas * means,
		var = betas * covs %*% betas
		)
 	
 )
}

