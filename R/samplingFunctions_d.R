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




#' Create a population with a specified mean and standard deviation
#' @param pop.mean Desired population mean
#' @param pop.sd d Desired population standard deviation
#' @param N Desired number of cases in population
#' @return Vector values defining the population
#' @examples
#' pop.1 <- make_population(pop.mean=100, pop.sd=15)
#' @export
make_population <- function(pop.mean=100, pop.sd=15,N=10000) {
     x <- as.numeric(scale(rnorm(N)))*pop.sd+pop.mean
     return(x)
}


#' Calculate a second population the differs by d from an input population
#' @param pop.data A vector of values defining population 1
#' @param d The amount by which population 2 should differ from population 1
#' @return Vector values for population 2
#' @examples
#' # Create two populations that differ by d = 1.0
#' pop.1 <- make_population(pop.mean=100, pop.sd=15)
#' pop.2 <- make_population_offset_by_d(pop.data=pop.1,d=1.0)
#' @export
make_population_offset_by_d <- function(pop.data,d) {
     n <- length(pop.data)
     x.m <- mean(pop.data)
     x.sd <- sd(pop.data)
     new.data <- (as.numeric(scale(rnorm(n))) + d)*x.sd + x.m
     return(new.data)
}



#' Calculate a number of sample d-values (unbiased) based on two populations
#' @param pop1 A vector of values defining population 1
#' @param pop2 A vector of values defining population 2
#' @param n Cell size for both cells for all samples. If you use n, do not use n.min or n.max.
#' @param n.min Minimum cell size across samples. Cell sizes will be equal in a single sample.
#' @param n.max Maximum cell size across samples. Cell sizes will be equal in a single sample.
#' @param number.of.samples Number of samples to obtain
#' @param number.of.decimals Number of decimals to report in returned data frame
#' @return Data frame with sample d-values
#' @examples
#' pop.1 <- make_population(pop.mean=100, pop.sd=15)
#' pop.2 <- make_population_offset_by_d(pop.data=pop.1,d=1.0)
#' my.samples <- get_d_samples_from_pops(pop1=pop.1,pop2=pop.2,n.min=50,n.max=150,number.of.samples=15)
#' @export
get_d_samples_from_pops <- function(pop1=NA,pop2=NA,n=NA,n.min=NA,n.max=NA,number.of.samples=10,number.of.decimals=2) {
     n.type <- "none"
     if (is.na(n)) {
          if (is.na(n.min)| is.na(n.max)){return()}
          #n.type <- "multiple"
          ns <- rep(n.min,number.of.samples) + round(runif(number.of.samples)*(n.max-n.min))
     } else {
          #n.type <- "single"
          ns <- rep(n,number.of.samples)
     }

     pop1.N <- length(pop1)
     pop2.N <- length(pop2)
     ds <- rep(NA,number.of.samples)
     for (i in 1:number.of.samples) {
          cur.n <- ns[i]
          id1 <- sample.int(pop1.N,cur.n)
          id2 <- sample.int(pop2.N,cur.n)
          group1.data <- pop1[id1]
          group2.data <- pop2[id2]
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
