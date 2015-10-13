#' Calculate x,y data to create a specified correlation
#' @param r Desired correlation
#' @param n Sample size
#' @return Data frame with sample data
#' @examples
#' my.data <- get_cor_data(r=.25, n=100)
#' @export
get_cor_data <- function(r, n) {
     x <- scale(matrix(rnorm(n*2),n,2),center=T,scale=T)
     svd.out <- svd(x)
     u <- svd.out$u

     y=u[,1]
     e=u[,2]

     if (round(r,3)!=0) {
          r2 <- r*r
          e.weight<-sqrt((1-r2)/r2)
          x = y + e.weight*e

          vec.out=matrix(NA,nrow=n,ncol=2)
          vec.out[,1]=scale(x,center=T,scale=T)
          vec.out[,2]=scale(y,center=T,scale=T)
     } else {
          vec.out=matrix(NA,nrow=n,ncol=2)
          vec.out[,1]=scale(y,center=T,scale=T)
          vec.out[,2]=scale(e,center=T,scale=T)
     }
     df.out <- data.frame(vec.out)
     names(df.out) <- c("x","y")
     return(df.out)
}











#' Calculate a number of sample correlations based on a specified population correlation
#' @param rho Population correlation. Do not use pop.data is you provide this value.
#' @param pop.data Population-level data frame for correlation. Do no use rho arguement if you provide this data.
#' @param n Sample size for all samples If you use n, do not use n.min or n.max.
#' @param n.min Minimum sample size
#' @param n.max Maximum sample size
#' @param number.of.samples Number of samples to obtain
#' @param number.of.decimals Number of decimals to report in returned data frame
#' @return Data frame with sample correlations
#' @examples
#' get_cor_samples(rho=.35,n=100)
#' my.samples <- get_cor_samples(rho=.35,n=100)
#' my.samples <- get_cor_samples(rho=.35,n.min=50,n.max=150,number.of.samples=15)
#' @export
get_cor_samples <- function(rho=NA,pop.data=NA,n=NA,n.min=NA,n.max=NA,number.of.samples=10,number.of.decimals=2) {

     n.type <- "none"
     if (is.na(n)) {
          if (is.na(n.min)| is.na(n.max)){return()}
          #n.type <- "multiple"
          ns <- rep(n.min,number.of.samples) + round(runif(number.of.samples)*(n.max-n.min))
     } else {
          #n.type <- "single"
          ns <- rep(n,number.of.samples)
     }

     pop.N <- 1000000 # Use a population with 1,000,000 people in it
     if (!is.na(rho)) {
          pop.data <- get_cor_data(r=rho,n=pop.N)
     } else {
          pop.N <- dim(pop.data)[1]
     }



     names(pop.data) <- c("x","y")
     rs <- rep(NA,number.of.samples)
     for (i in 1:number.of.samples) {
          cur.n <- ns[i]
          ids <- sample.int(pop.N,cur.n)
          sample.data <- pop.data[ids,]
          rs[i] <- round(stats::cor(x=sample.data$x,y=sample.data$y),number.of.decimals)
     }
     n <- ns
     r <- rs
     xx<-1:number.of.samples
     study.number <- paste("Study",xx)
     data.out <- data.frame(study.number,n,r)
     rownames(data.out) <- NULL

     return(data.out)
}



#' Sort data frame based on column r
#' @param df.in Data frame to sort
#' @return Data frame that has been sorted by column r
#' @examples
#' my.samples <- get_cor_samples(rho=.35,n.min=50,n.max=150,number.of.samples=10)
#' my.samples.sorted <- sort_samples_by_r(my.samples)
#' @export
sort_samples_by_r <- function(df.in){
     df.in <- df.in[order(df.in$r),]
     rownames(df.in) <-NULL
     return(df.in)
}


#' Calculate a data frame where columns correlate as specified by a correlation matrix
#' @param r.matrix Desired correlation matrix
#' @param n Sample size
#' @return Data frame with sample data
#' @examples
#' N <- 600
#' M <- matrix(c(1.00, 0.60, 0.30, 0.30,
#'             0.60, 1.00, 0.00, 0.60,
#'             0.30, 0.00, 1.00, 0.00,
#'             0.30, 0.60, 0.00, 1.00), nrow=4, ncol=4)
#' my.data <- get_cor_data_from_matrix(M,N)
#' library(apaTables)
#' apa.cor.table(my.data)
#' @export
get_cor_data_from_matrix <- function(r.matrix,n) {
     L <-  chol(r.matrix)
     tL <-t(L)
     nvars <- dim(L)[1]

     x=scale(matrix(rnorm(n*nvars),n,nvars),center=T,scale=T)
     svd.out <- svd(x)
     u <- svd.out$u

     matrix.out <- u %*% L
     df.out <- data.frame(matrix.out)
     return(df.out)
}

#' Change the mean and standard deviation of each column in a data frame to specified values
#' @param cor.data Data frame
#' @param means.in A vector of desired column means. The number means must correspond to the number of columns.
#' @param sds.in A vector of desired column standard deviations. The number standard deviations must correspond to the number of columns.
#' @return Data frame with the desired means and standard deviations for the columns
#' @examples
#' N <- 600
#' M <- matrix(c(1.00, 0.60, 0.30, 0.30,
#'             0.60, 1.00, 0.00, 0.60,
#'             0.30, 0.00, 1.00, 0.00,
#'             0.30, 0.60, 0.00, 1.00), nrow=4, ncol=4)
#' my.data <- get_cor_data_from_matrix(M,N)
#' desired.means <- c(3.5,4,4.5,5)
#' desired.sds <- c(1.1,1.2,1.3,1.4)
#' my.data <- cor_data_set_m_sd(my.data, desired.means,desired.sds)
#' library(apaTables)
#' apa.cor.table(my.data)
#' @export
cor_data_set_m_sd <- function(cor.data,means.in,sds.in) {
     df.out <- cor.data
     n.col <- ncol(df.out)

     for (cur.col.num in 1:n.col) {
          cur.col <- df.out[,cur.col.num]
          cur.col <- as.numeric(scale(cur.col))
          cur.col <- cur.col*sds.in[cur.col.num]+means.in[cur.col.num]
          df.out[,cur.col.num]<-cur.col
     }
     return(round(df.out,2))
}

