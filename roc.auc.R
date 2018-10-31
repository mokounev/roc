roc.auc <- function(df, ess.genes){
  
  library(zoo)
  
  ## Check input
  if (!is.data.frame(df)) {
    stop("Please provide a data.frame for the df arguement")
  }
  if (!"gene" %in% colnames(df)) {
    stop("Please provide df with column 'gene'")
  }
  if (!"rho.depletion.gene" %in% colnames(df)){3
    stop("Please provide df with column 'rho.depletion.gene'")
  }
  if (!is.character(ess.genes)) {
    stop("Variable 'essential.genes' must be of type 'character'")
  } 
  
  ## Remove duplicate genes and order df by rho depletion score
  df <- df[!duplicated(df$gene), ]
  df <- df[order(df$rho.depletion.gene), ]
  
  label <- as.integer(df$gene %in% ess.genes)
  score <- seq_len(nrow(df))

  label <- label[order(score, decreasing = TRUE)]  

  result <- data.frame(TPR = 1 - (cumsum(label)/sum(label)), FPR = 1 - (cumsum(!label)/sum(!label)), label)

  x <- result$TPR
  y <- result$FPR
  
  if (length(x) != length(y))
    stop("length x must be equal length y")
  y <- y[x >= 0 & x <= 1]
  x <- x[x >= 0 & x <= 1]
  id <- order(x)
  auc <-  sum(diff(x[id])*rollmean(y[id],2))
  
  plot(x, y, type = 'l', xlab = "TPR", ylab = "FPR")
  auc

}

