
#' @title exp2lat
#' @description Transforms means and covariance matrix from expected (log-normal) to latent (normal) scale
#' @param mean a vector of means. These need to all be >0.
#' @param cov a single variance or a covariance matrix
#' @details Transforms expected mean(s) and (co)variance(s) of a multivariate normal distribution to expected mean(s) and (co)variance(s) of a multivariate log-normal distribution. This is useful when trying to go from expected (log-normal) to expected (log-normal) scales in a Poisson GLM with a log link function.
#' @author Joel Pick - joel.l.pick@gmail.com
#' @return a list, which an element 'mean' for the transformed mean(s) and another element 'cov' for the transformed covariance matrix
#' @examples
#' exp2lat(1,1)
#' 
#' @export
#' 
exp2lat<- function(mean,cov){

	if(!(is.numeric(mean) & is.numeric(mean)))  stop( "Inputs must be numeric", call.=FALSE)
	
	if(!is.vector(mean))  stop( "mean must be vector", call.=FALSE)
	if(any(mean<=0)) stop("elements of mean should all be >0", call.=FALSE)
	mean <- as.matrix(mean) #make a column vector

	if(!(is.matrix(cov) | is.vector(cov)))  stop( "cov must be a vector or matrix", call.=FALSE)
	if(is.matrix(cov)){
		if(ncol(cov)!=nrow(cov)) stop('cov must be a square matrix', call.=FALSE)
		if(!isSymmetric(cov)) stop("cov must be a symmetric matrix", call.=FALSE)
		if(length(mean)!=nrow(cov)) stop('mean must be the same length as the number of rows/columns in cov', call.=FALSE)
		if(any(eigen(cov)$values<0))stop("cov matrix should be positive semi-definite", call.=FALSE)	
	}
	if(is.vector(cov)){
		if(length(mean)!=length(cov)) stop('mean must be the same length as cov', call.=FALSE)
		cov <- diag(cov, nrow=length(cov),ncol=length(cov)) 
	}
	
	mean_out <- rep(NA,nrow(cov))
	for(i in 1:nrow(cov)) mean_out[i] <- log(mean[i]^2/sqrt(mean[i]^2+cov[i,i]))
	
	cov_out <- matrix(NA,nrow(cov),ncol(cov))
	for(i in 1:nrow(cov)){
		for(j in 1:ncol(cov)) cov_out[i,j] <- log(1 + cov[i,j]/(mean[i]*mean[j]))
	}
	return(list(mean=mean_out,cov=cov_out))
}


#' @title lat2exp
#' @description Transforms means and covariance matrix from latent (normal) to expected (log-normal) scale
#' @param mean a vector of means
#' @param cov a single variance or a covariance matrix
#' @details Transforms expected mean(s) and (co)variance(s) of a multivariate log-normal distribution to expected mean(s) and (co)variance(s) of a multivariate normal distribution. This is useful when trying to go from latent (normal) to expected (log-normal) scales in a Poisson GLM with a log link function.
#' @author Joel Pick - joel.l.pick@gmail.com
#' @return a list, which an element 'mean' for the transformed mean(s) and another element 'cov' for the transformed covariance matrix
#' @examples
#' lat2exp(1,1)
#' 
#' @export

lat2exp<- function(mean,cov){
	
	if(!(is.numeric(mean) & is.numeric(mean)))  stop( "Inputs must be numeric", call.=FALSE)
	
	if(!is.vector(mean))  stop( "mean must be vector", call.=FALSE)
	mean <- as.matrix(mean) #make a column vector

	if(!(is.matrix(cov) | is.vector(cov)))  stop( "cov must be a vector or matrix", call.=FALSE)
	if(is.matrix(cov)){
		if(ncol(cov)!=nrow(cov)) stop('cov must be a square matrix', call.=FALSE)
		if(!isSymmetric(cov)) stop("cov must be a symmetric matrix", call.=FALSE)
		if(length(mean)!=nrow(cov)) stop('mean must be the same length as the number of rows/columns in cov', call.=FALSE)
		if(any(eigen(cov)$values<0))stop("cov matrix should be positive semi-definite", call.=FALSE)	
	}
	if(is.vector(cov)){
		if(length(mean)!=length(cov)) stop('mean must be the same length as cov', call.=FALSE)
		cov <- diag(cov, nrow=length(cov),ncol=length(cov)) 
	}
	
	mean_out <- rep(NA,nrow(cov))
	for(i in 1:nrow(cov)) mean_out[i] <- exp(mean[i]+ cov[i,i]/2)

	cov_out <- matrix(NA,nrow(cov),ncol(cov))
	for(i in 1:nrow(cov)){
		for(j in 1:ncol(cov)) cov_out[i,j] <- exp(mean[i]+mean[j] + (cov[i,i]+cov[j,j])/2)*(exp(cov[i,j])-1)
	}
	return(list(mean=mean_out,cov=cov_out))
}