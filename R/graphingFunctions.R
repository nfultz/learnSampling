
create_numbers_within_bins <- function(df, center_bins) {
     bin_number <- df$bin_number
     max_bin_number <- max(bin_number)
     position_in_bin <- rep(NA, length(bin_number))
     bin_center <- rep(NA, length(bin_number))
     for (cur_bin in 1:max_bin_number) {
          is_cur_bin  <- cur_bin==bin_number
          num_cur_bin <- sum(is_cur_bin)
          if (num_cur_bin>0) {
               position_in_cur_bin <- 1:num_cur_bin
               position_in_bin[is_cur_bin] <- position_in_cur_bin
               bin_center[is_cur_bin] <- rep(center_bins[cur_bin], num_cur_bin)
          }
     }

     df$position_in_bin <- position_in_bin
     df$bin_center  <- bin_center
     return(df)
}

make_graph <- function(df, pop.value, plot.type) {
     x <- sort(df$x)
     N <- length(x)

     if (plot.type=="r") {
          pop_str <- paste("rho == ", pop.value)
          plot_char <- rep("r",N)

     } else if (plot.type=="d") {
          pop_str <- paste("delta == ", pop.value)
          plot_char <- rep("d",N)

     } else if (plot.type=="M") {
          pop_str <- paste("mu == ", pop.value)
          plot_char <- rep("M",N)
     } else {
          pop_str <- paste("mu == ", pop.value)
          plot_char <- rep("X",N)
     }



     number_bins <- 11
     x_min <- min(x)
     x_max <- max(x)
     x_mean <- mean(x)

     x_diff     <- max(c((x_mean-x_min),(x_max-x_mean)))
     bin_width  <- x_diff/(number_bins/2)
     bin_min    <- x_mean - (number_bins/2)*bin_width
     bin_max    <- x_mean + (number_bins/2)*bin_width
     bin_breaks <- seq(bin_min,bin_max, by=bin_width)

     lower_bin_bounds <- bin_breaks[1:(length(bin_breaks)-1)]
     upper_bin_bounds <- bin_breaks[2:(length(bin_breaks))]
     temp_bin_bounds  <- rbind(lower_bin_bounds, upper_bin_bounds)
     center_bins      <- colMeans(temp_bin_bounds)

     bin_number    <- cut(x, breaks=bin_breaks,labels=FALSE)
     bin_number[1] <- 1

     dout <- df
     dout <- cbind(df, bin_number)
     dout <- create_numbers_within_bins(dout, center_bins)
     dout$position_in_bin <- dout$position_in_bin

     y_max_value <- max(dout$position_in_bin) * 1.20
     y_axis_max <- y_max_value * 1.30

     dout$pchar <- plot_char

     library(ggplot2)
     dplot <- ggplot2::ggplot(data=dout, aes(x=bin_center, y=position_in_bin, label=pchar))
     dplot <- dplot + ggplot2::coord_cartesian(ylim=c(0,y_axis_max))
     dplot <- dplot + ggplot2::geom_text(check_overlap = TRUE, lineheight=1.5, fontface="bold.italic", hjust=0.5)
     dplot <- dplot + ggplot2::geom_segment(aes(x=pop.value, y=0,xend=pop.value,yend=y_max_value),size=.1)
     dplot <- dplot + ggplot2::annotate("text", x = pop.value, y = y_max_value*1.05, label = pop_str, parse=TRUE)
     dplot <- dplot + ggplot2::labs(x="Sample Value", y="Frequency")
     dplot <- dplot + ggplot2::theme_classic()

     return(dplot)
}


#df <- data.frame(x=rnorm(500, mean = 100, sd=15))
#gg <- make_graph(df, pop.value = 100, plot.type="")
#print(gg)
