

#cohensD ---------------------------------
cohend_gw <- function(x, y, hedges) {
  
  cohend <- function(x, y) {
    cohend_est = try(cohen.d(
      as.numeric(x),
      as.numeric(y),
      na.rm = TRUE,
      pooled = FALSE,
      paired = FALSE,
      hedges.correction = hedges
    ))
    if (is(cohend_est, "try-error"))
      return(NA)
    else
      return(signif((cohend_est$estimate), digits = 3))
  }
  
  x <- data.frame(t(x))
  y <- data.frame(t(y))
  
  cd <- mapply(function(x,y) cohend(x,y), x, y)
  return(cd)
}

#missing factor ---------------------------------
missing_factor_gw <- function(x, y) {
  x <- x %>% mutate_all(as.numeric)
  y <- y %>% mutate_all(as.numeric)
  mf_x <- rowSums(x) / ncol(x)
  mf_y <- rowSums(y) / ncol(y)
  df_mf <- data.frame(cbind(mf_x, mf_y), stringsAsFactors = FALSE)
  df_mf$max <-
    apply(
      df_mf,
      1,
      FUN = function(x) {
        max(x, na.rm = TRUE)
      }
    )
  return(signif(df_mf$max, digits = 3))
}


#fold change ---------------------------------
foldchange_gw <- function(x, y) {
  if (!as.logical(dpmsr_set$x$pair_comp)) {
    ave_x = rowMeans(x)
    ave_y = rowMeans(y)
    test = ave_x / ave_y
  } else{
    sn <- ncol(x)
    indiv_fc <- x
    for (i in 1:sn) {
      indiv_fc[i] <- (x[i] / y[i])
    }
    test <- rowMeans(indiv_fc)
  }
  fc <- ifelse ((test >= 1), test,-1 / test)
  return(signif(fc, digits = 7))
}

#fold change pair---------------------------------
foldchange_pair_gw <- function(x, y) {
  sn <- ncol(x)
  indiv_fc <- x
  for (i in 1:sn) {
    indiv_fc[i] <- (x[i] / y[i])
  }
  test <- rowMeans(indiv_fc)
  fc <- ifelse ((test >= 1), test,-1 / test)
  return(signif(fc, digits = 7))
}



#fold change decimal ---------------------------------
foldchange_decimal_gw <- function(x, y) {
  if (!as.logical(dpmsr_set$x$pair_comp)) {
    ave_x = rowMeans(x)
    ave_y = rowMeans(y)
    test = ave_x / ave_y
  } else{
    sn <- ncol(x)
    indiv_fc <- x
    for (i in 1:sn) {
      indiv_fc[i] <- (x[i] / y[i])
    }
    test <- rowMeans(indiv_fc)
  }
  fc <- test
  return(signif(fc, digits = 7))
}

#fold change pair decimal---------------------------------
foldchange_pair_decimal_gw <- function(x, y) {
  sn <- ncol(x)
  indiv_fc <- x
  for (i in 1:sn) {
    indiv_fc[i] <- (x[i] / y[i])
  }
  test <- rowMeans(indiv_fc)
  fc <- test
  return(signif(fc, digits = 7))
}


#Percent CV ---------------------------------
percentCV_gw <- function(x) {
  ave <- rowMeans(x)
  n <- ncol(x)
  sd <- apply(x[1:n], 1, sd)
  cv <- (100 * sd / ave)
  return(signif(cv, digits = 3))
}


#t.test ---------------------------------
ttest_gw <- function(x, y) {
  if (!as.logical(dpmsr_set$x$pair_comp)) {
    ttest_pvalue = try(t.test(
      x,
      y,
      alternative = "two.sided",
      var.equal = FALSE,
      paired = FALSE
    ),
    silent = TRUE)
  } else{
    ttest_pvalue = t.test(x, y, paired = TRUE)
  }
  if (is(ttest_pvalue, "try-error"))
    return(NA)
  else
    return(signif((ttest_pvalue$p.value), digits = 7))
}


pvalue_gw <- function(x, y) {
  x <- log2(x)
  y <- log2(y)
  temp_pval <- rep(NA, nrow(x))
  for (i in 1:nrow(x))
  {
    temp_pval[i] <- ttest_gw(as.numeric(x[i, ]), as.numeric(y[i, ]))
  }
  return(temp_pval)
}



exactTest_gw <- function(x, y) {
  #x <- log2(x)
  #y <- log2(y)
  et <- exactTestDoubleTail(x, y)
  return(et)
}



#x<-comp_N_data
#y<-comp_D_data
#comp_name <- comp_groups$comp_name[1]
limma_gw <- function(x, y) {   #, comp_name, plot_dir) {
  xy <- cbind(x, y)
  xy <- log2(xy)
  n <- ncol(x)
  d <- ncol(y)
  design <- model.matrix( ~ 0 + factor(c(rep(1, n), rep(0, d))))
  colnames(design) <- c("group1", "group2")
  contrast.matrix <- makeContrasts(group2 - group1, levels = design)
  fit <- lmFit(xy, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  topfit <- topTable(fit2,
                     coef = 1,
                     sort = "none",
                     number = Inf)
  data_out <- topfit$P.Value
  #try(limma_qq(fit2$t, comp_name, plot_dir), silent = TRUE)
  #try(limma_ma(topfit, comp_name, plot_dir), silent = TRUE)
  #try(limma_volcano(fit2, comp_name, plot_dir), silent = TRUE)
  return(data_out)
}


old_limma_gw <- function(x, y) {
  xy <- cbind(x, y)
  xy <- log2(xy)
  n <- ncol(x)
  d <- ncol(y)
  design <- cbind(Grp1 = 1, Grp2vs1 = c(rep(1, n), rep(0, d)))
  #design <- c(0,0,0,1,1,1)
  # Ordinary fit
  fit <- lmFit(xy, design)
  fit <- eBayes(fit)
  topfit <-
    topTable(
      fit,
      coef = 2,
      adjust.method = "BH",
      sort = "none",
      number = Inf
    )
  data_out <- topfit$P.Value
  return(data_out)
}
