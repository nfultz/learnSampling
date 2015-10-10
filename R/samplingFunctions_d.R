#' Calculate a number of sample d-values (unbiased) based on a specified population correlation
#' @param pop.d Population d-value
#' @param n Cell size for both cells for all samples. If you use n, do not use n.min or n.max.
#' @param n.min Minimum cell size across samples. Cell sizes will be equal in a single sample.
#' @param n.max Maximum cell size across samples. Cell sizes will be equal in a single sample.
#' @param number.of.samples Number of samples to obtain
#' @param number.of.decimals Number of decimals to report in returned data frame
#' @return Data frame with sample d-values
#' @examples
#' get_d_samples(pop.d=.35,n=100)
#' my.samples <- get_d_samples(pop.d=.35,n=100)
#' my.samples <- get_d_samples(pop.d=.35,n.min=50,n.max=150,number.of.samples=15)
#' @export
get_d_samples <- function(pop.d=NA,n=NA,n.min=NA,n.max=NA,number.of.samples=10,number.of.decimals=2) {

     if (is.na(pop.d)) {return()}

     n.type <- "none"
     if (is.na(n)) {
          if (is.na(n.min)| is.na(n.max)){return()}
          #n.type <- "multiple"
          ns <- rep(n.min,number.of.samples) + round(runif(number.of.samples)*(n.max-n.min))
     } else {
          #n.type <- "single"
          ns <- rep(n,number.of.samples)
     }

     ds <- rep(NA,number.of.samples)
     for (i in 1:number.of.samples) {
          cur.n <- ns[i]
          group1.data <- rnorm(cur.n) + pop.d
          group2.data <- rnorm(cur.n)
          cur.d <- MBESS::smd(Group.1=group1.data,Group.2=group2.data,Unbiased=TRUE)
          ds[i] <- round(cur.d, number.of.decimals)
     }
     n <- ns
     d <- ds
     xx<-1:number.of.samples
     study.number <- paste("Study",xx)
     data.out <- data.frame(study.number,n,d)
     rownames(data.out) <- NULL

     return(data.out)
}

#' Sort data frame based on column d
#' @param df.in Data frame to sort
#' @return Data frame that has been sorted by column d
#' @examples
#' my.samples <- get_d_samples(pop.d=.8,n.min=50,n.max=150,number.of.samples=10)
#' my.samples.sorted <- sort_samples_by_d(my.samples)
#' @export
sort_samples_by_d <- function(df.in){
     df.in <- df.in[order(df.in$d),]
     rownames(df.in) <-NULL
     return(df.in)
}

