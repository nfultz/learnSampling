#' Calculate x,y data to create a specified correlation
#' @param r Desired correlation
#' @param n Sample size
#' @return Matrix with sample data
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
     return(df.out)
}

#' Calculate a number of sample correlations based on a specified population correlation
#' @param rho Population correlation
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
get_cor_samples <- function(rho=NA,n=NA,n.min=NA,n.max=NA,number.of.samples=10,number.of.decimals=2) {

     if (is.na(rho)) {return()}

     n.type <- "none"
     if (is.na(n)) {
          if (is.na(n.min)| is.na(n.max)){return()}
          #n.type <- "multiple"
          ns <- rep(n.min,number.of.samples) + round(runif(number.of.samples)*(n.max-n.min))
     } else {
          #n.type <- "single"
          ns <- rep(n,number.of.samples)
     }

     pop.N <- 1000000
     pop.data <- get_cor_data(r=rho,n=pop.N)
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
