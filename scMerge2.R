
library(rsvd)
library(BiocSingular)

# source("/albona/nobackup2/yingxinl/scMerge2/scMerge2_functions/construct_pseudoBulk.R")

# 
# residop_fast <-
#   function (A, B)
#   {
#     return(A - B %*% solve(t(B) %*% B) %*% (t(B) %*% A))
#   }

# tological <-
#   function (ctl, n) 
#   {
#     ctl2 = rep(FALSE, n)
#     ctl2[ctl] = TRUE
#     return(ctl2)
#   }

tological <- function(ctl, n) {
  ctl2 <- rep(FALSE, n)
  ctl2[ctl] <- TRUE
  return(ctl2)
}


#' @importFrom DelayedArray t
my_residop <- function(A, B){
  tBB = DelayedArray::t(B) %*% B
  tBB_inv = Matrix::solve(tBB)
  BtBB_inv = B %*% tBB_inv
  tBA = DelayedArray::t(B) %*% A
  BtBB_inv_tBA = BtBB_inv %*% tBA
  return(A - BtBB_inv_tBA)
}





####################################################### 
zeroOneScale <- function(v) {
  v <- (v + 1)/2
  return(v)
}
####################################################### 
standardize <- function(exprsMat, batch) {
  num_cell <- ncol(exprsMat)
  num_batch <- length(unique(batch))
  batch <- as.factor(batch)
  grand.mean <- matrix(base::rowMeans(exprsMat), nrow = 1)
  stand.mean <- t(grand.mean) %*% t(rep(1, num_cell))
  design <- stats::model.matrix(~-1 + batch)
  B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(exprsMat)))
  var.pooled <- ((exprsMat - t(design %*% B.hat))^2) %*% rep(1/(num_cell - 
                                                                  num_batch), num_cell)
  s.data <- (exprsMat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, 
                                                                num_cell)))
  return(res = list(s.data = s.data, stand.mean = stand.mean, 
                    stand.var = var.pooled))
}
###############################
standardize2 <- function(Y, batch) {
  num_cell <- ncol(Y)
  num_batch <- length(unique(batch))
  batch <- as.factor(batch)
  stand_mean <- DelayedArray::rowMeans(Y)
  design <- stats::model.matrix(~-1 + batch)
  B.hat = solve_axb(a = DelayedArray::t(design) %*% design,
                    b = DelayedArray::t(Y %*% design))
  
  stand_var <- DelayedArray::rowSums(((Y - DelayedArray::t(B.hat) %*% DelayedArray::t(design))^2))/(num_cell - num_batch)
  stand_Y <- (Y-stand_mean)/sqrt(stand_var)
  return(res = list(stand_Y = stand_Y, 
                    stand_mean = stand_mean, 
                    stand_var = stand_var))
}
####################################################### 
f_measure <- function(cell_type, batch) {
  f <- 2 * (cell_type * batch)/(cell_type + batch)
  return(f)
}
####################################################### 
calculateSil <- function(x, fast_svd, cell_type, batch) {
  if (fast_svd & !any(dim(x$newY) < 50)) {
    pca.data <- irlba::prcomp_irlba(x$newY, n = 10)
  } else {
    pca.data <- stats::prcomp(x$newY)
  }
  
  result = c(kBET_batch_sil(pca.data, as.numeric(as.factor(cell_type)), 
                            nPCs = 10), kBET_batch_sil(pca.data, as.numeric(as.factor(batch)), 
                                                       nPCs = 10))
  return(result)
}

kBET_batch_sil <- function(pca.data, batch, nPCs = 10) {
  ## This function was copied from kBET, which cannot be
  ## imported because it is only a GitHub package
  dd <- as.matrix(stats::dist(pca.data$x[, seq_len(nPCs)]))
  score_sil <- summary(cluster::silhouette(as.numeric(batch), 
                                           dd))$avg.width
  return(score_sil)
}




###############################
solve_axb = function(a, b){
  x = solve(DelayedArray::t(a) %*% a) %*% DelayedArray::t(a) %*% b
  return(x)
}

fastRUVIII <-  function(Y, M, ctl, k = NULL, eta = NULL, 
                        svd_k = 50,
                        include.intercept = TRUE, batch = NULL,
                        average = FALSE, fullalpha = NULL, return.info = FALSE, inputcheck = TRUE,
                        regularise = 0.01, scale = TRUE, 
                        BPPARAM = SerialParam(), BSPARAM = ExactParam()) 
{
  if (scale) {
    scale_res <- standardize2(t(Y), batch)
    stand_tY <- DelayedArray::t(scale_res$stand_Y)
    stand_sd <- sqrt(scale_res$stand_var)
    stand_mean <- scale_res$stand_mean
    Y <- stand_tY
    
    rm("stand_tY")
  } 
  
  if (is.data.frame(Y)) 
    Y = data.matrix(Y)
  m = nrow(Y)
  n = ncol(Y)
  M = replicate.matrix(M)
  ctl = tological(ctl, n)
  if (inputcheck) {
    if (m > n) 
      warning("m is greater than n!  This is not a problem itself, but may indicate that you need to transpose your data matrix.  Please ensure that rows correspond to observations (e.g. microarrays) and columns correspond to features (e.g. genes).")
    if (sum(is.na(Y)) > 0) 
      warning("Y contains missing values.  This is not supported.")
    if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 
        0) 
      warning("Y contains infinities.  This is not supported.")
  }
  Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
  
  
  if (class(BSPARAM) != "ExactParam") {
    svd_k <- min(m - ncol(M), sum(ctl), svd_k, na.rm = TRUE)
  } else {
    svd_k <- min(m - ncol(M), sum(ctl), na.rm = TRUE)
  }
  
  ## m represent the number of samples/observations ncol(M)
  ## represent the number of replicates If the replicate matrix
  ## is such that we have more replicates than samples, then
  ## RUV3 is not appropriate, thus, we return the Original input
  ## matrix
  if (ncol(M) >= m | k == 0) {
    newY <- Y
    fullalpha <- NULL
  } else {
    
    if (is.null(fullalpha)) 
    {
      ## The main RUVIII process Applies the residual operator of a
      ## matrix M to a matrix Y Y0 has the same dimensions as Y,
      ## i.e. m rows (observations) and n columns (genes).
      
      # if(class(Y) == "matrix"){
      #   Y0 <- eigenResidop(Y, M)
      # } else if (class(Y) == "dgeMatrix"){
      #   Y0 <- eigenResidop(as.matrix(Y), M)
      # } else {
      #   Y0 <- my_residop(Y, M)
      # }
      
      Y0 <- my_residop(Y, M)
      
      svdObj <- BiocSingular::runSVD(x = Y0, k = svd_k, BPPARAM = BPPARAM, BSPARAM = BSPARAM)
      
      fullalpha <- DelayedArray::t(svdObj$u[, seq_len(svd_k), drop = FALSE]) %*% Y
    }  ## End is.null(fullalpha)
    ###############
    alpha <- fullalpha[seq_len(min(k, nrow(fullalpha))), , drop = FALSE]
    ac <- alpha[, ctl, drop = FALSE]
    W <- Y[, ctl] %*% DelayedArray::t(ac) %*% solve(ac %*% DelayedArray::t(ac))
    newY <- Y - W %*% alpha
  }  ## End else(ncol(M) >= m | k == 0)
  
  
  
  
  if (!return.info) 
    return(newY)
  else return(list(newY = newY, M = M, fullalpha = fullalpha, W = W))
}





Pm <- function(B) {
  tBB = DelayedArray::t(B) %*% B
  tBB_inv = Matrix::solve(tBB)
  BtBB_inv = B %*% tBB_inv
  res = BtBB_inv %*% DelayedArray::t(B)
  return(res)
}

residop2 <- function(Pm, Y) {
  res = Y - Pm %*% Y
  return(res)
}


fastRUVIII2 <-  function(Y, Pm, ctl, k = NULL, eta = NULL, 
                         svd_k = 50,
                         include.intercept = TRUE, batch = NULL,
                         average = FALSE, fullalpha = NULL, return.info = FALSE, inputcheck = TRUE,
                         regularise = 0.01, scale = TRUE, 
                         BPPARAM = SerialParam(), BSPARAM = ExactParam()) 
{
  if (scale) {
    scale_res <- standardize2(t(Y), batch)
    stand_tY <- DelayedArray::t(scale_res$stand_Y)
    stand_sd <- sqrt(scale_res$stand_var)
    stand_mean <- scale_res$stand_mean
    Y <- stand_tY
    
    rm("stand_tY")
  } 
  
  if (is.data.frame(Y)) 
    Y = data.matrix(Y)
  m = nrow(Y)
  n = ncol(Y)
  # M = replicate.matrix(M)
  ctl = tological(ctl, n)
  if (inputcheck) {
    if (m > n) 
      warning("m is greater than n!  This is not a problem itself, but may indicate that you need to transpose your data matrix.  Please ensure that rows correspond to observations (e.g. microarrays) and columns correspond to features (e.g. genes).")
    if (sum(is.na(Y)) > 0) 
      warning("Y contains missing values.  This is not supported.")
    if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 
        0) 
      warning("Y contains infinities.  This is not supported.")
  }
  Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
  
  
  if (class(BSPARAM) != "ExactParam") {
    #svd_k <- min(m - ncol(M), sum(ctl), svd_k, na.rm = TRUE)
    svd_k <- min(m, sum(ctl), svd_k, na.rm = TRUE)
  } else {
    #svd_k <- min(m - ncol(M), sum(ctl), na.rm = TRUE)
    svd_k <- min(m, sum(ctl), na.rm = TRUE)
  }
  
  ## m represent the number of samples/observations ncol(M)
  ## represent the number of replicates If the replicate matrix
  ## is such that we have more replicates than samples, then
  ## RUV3 is not appropriate, thus, we return the Original input
  ## matrix
  if (is.null(fullalpha)) 
  {
    
    
    #Y0 <- my_residop(Y, M)
    Y0 <- residop2(Pm, Y)
    
    svdObj <- BiocSingular::runSVD(x = Y0, k = svd_k, BPPARAM = BPPARAM, BSPARAM = BSPARAM)
    
    fullalpha <- DelayedArray::t(svdObj$u[, seq_len(svd_k), drop = FALSE]) %*% Y
  }  ## End is.null(fullalpha)
  ###############
  alpha <- fullalpha[seq_len(min(k, nrow(fullalpha))), , drop = FALSE]
  ac <- alpha[, ctl, drop = FALSE]
  W <- Y[, ctl] %*% DelayedArray::t(ac) %*% solve(ac %*% DelayedArray::t(ac))
  newY <- Y - W %*% alpha
  
  if (!return.info) 
    return(newY)
  else return(list(newY = newY, Pm = Pm, fullalpha = fullalpha, W = W))
}





# 
# 
# 
# pseudoRUVIII <-  function(Y, Y_pseudo, M, ctl, k = NULL, eta = NULL, include.intercept = TRUE, batch = NULL,
#                           average = FALSE, fullalpha = NULL, return.info = FALSE, inputcheck = TRUE,
#                           regularise = 0.01, scale = TRUE, 
#                           BPPARAM = SerialParam(), BSPARAM = ExactParam()) 
# {
#   # if (scale) {
#   #   scale_res <- standardize2(t(Y), batch)
#   #   stand_tY <- DelayedArray::t(scale_res$stand_Y)
#   #   stand_sd <- sqrt(scale_res$stand_var)
#   #   stand_mean <- scale_res$stand_mean
#   #   Y <- stand_tY
#   #   
#   #   rm("stand_tY")
#   # } 
#   
#   if (is.data.frame(Y)) 
#     Y = data.matrix(Y)
#   m = nrow(Y)
#   n = ncol(Y)
#   
#   if (is.data.frame(Y_pseudo)) 
#     Y = data.matrix(Y_pseudo)
#   m_pseudo = nrow(Y_pseudo)
#   n_pseudo = ncol(Y_pseudo)
#   
#   
#   M = replicate.matrix(M)
#   ctl = tological(ctl, n)
#   if (inputcheck) {
#     if (m > n) 
#       warning("m is greater than n!  This is not a problem itself, but may indicate that you need to transpose your data matrix.  Please ensure that rows correspond to observations (e.g. microarrays) and columns correspond to features (e.g. genes).")
#     if (sum(is.na(Y)) > 0) 
#       warning("Y contains missing values.  This is not supported.")
#     if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 
#         0) 
#       warning("Y contains infinities.  This is not supported.")
#   }
#   Y_pseudo = RUV1(Y_pseudo, eta, ctl, include.intercept = include.intercept)
#   
#   
#   if (class(BSPARAM) != "ExactParam") {
#     svd_k <- min(m_pseudo - ncol(M), sum(ctl), svd_k, na.rm = TRUE)
#   } else {
#     svd_k <- min(m_pseudo - ncol(M), sum(ctl), na.rm = TRUE)
#   }
#   
#   ## m represent the number of samples/observations ncol(M)
#   ## represent the number of replicates If the replicate matrix
#   ## is such that we have more replicates than samples, then
#   ## RUV3 is not appropriate, thus, we return the Original input
#   ## matrix
#   if (ncol(M) >= m_pseudo | k == 0) {
#     newY <- Y_pseudo
#     fullalpha <- NULL
#   } else {
#     
#     if (is.null(fullalpha)) 
#     {
#       ## The main RUVIII process Applies the residual operator of a
#       ## matrix M to a matrix Y Y0 has the same dimensions as Y,
#       ## i.e. m rows (observations) and n columns (genes).
#       
#       # if(class(Y) == "matrix"){
#       #   Y0 <- eigenResidop(Y_pseudo, M)
#       # } else if (class(Y) == "dgeMatrix"){
#       #   Y0 <- eigenResidop(as.matrix(Y_pseudo), M)
#       # } else {
#       #   Y0 <- my_residop(Y_pseudo, M)
#       # }
#       
#       Y0 <- my_residop(Y_pseudo, M)
#       svdObj <- BiocSingular::runSVD(x = Y0, k = svd_k, BPPARAM = BPPARAM, BSPARAM = BSPARAM)
#       
#       fullalpha <- DelayedArray::t(svdObj$u[, seq_len(svd_k), drop = FALSE]) %*% Y_pseudo
#     }  ## End is.null(fullalpha)
#     ###############
#     alpha <- fullalpha[seq_len(min(k, nrow(fullalpha))), , drop = FALSE]
#     ac <- alpha[, ctl, drop = FALSE]
#     W <- Y[, ctl] %*% DelayedArray::t(ac) %*% solve(ac %*% DelayedArray::t(ac))
#     print(dim(alpha))
#     print(dim(W))
#     
#     
#     if (nrow(Y) > 50000 & ncol(Y) > 10000) {
#       mid <- round(ncol(alpha)/2)
#       print(mid)
#       newY <- Y[, seq_len(mid)] - W %*% alpha[, seq_len(mid)]
#       newY2 <- Y[, c((mid + 1):ncol(alpha))] - W %*% alpha[, c((mid + 1):ncol(alpha))]
#       return(list(newY, newY2))
#     }  else {
#       newY <- Y - W %*% alpha
#     }
#     
#   }  ## End else(ncol(M) >= m | k == 0)
#   
#   rm("Y_pseudo")
#   
#   
#   if (!return.info) 
#     return(newY)
#   else return(list(newY = newY, M = M, fullalpha = fullalpha))
# }
# 



pseudoRUVIII <- function(Y, Y_pseudo, M, ctl, k = NULL, eta = NULL, 
                         include.intercept = TRUE, inputcheck = TRUE, 
                         #batch = NULL, regularise = 0.01, scale = TRUE, 
                         #average = FALSE, 
                         fullalpha = NULL, return.info = FALSE, 
                         subset = NULL,
                         BPPARAM = SerialParam(), 
                         BSPARAM = ExactParam(),
                         normalised = TRUE) 
{
  
  
  
  
  m = nrow(Y)
  n = ncol(Y)
  
  if (is.data.frame(Y_pseudo)) {
    Y = data.matrix(Y_pseudo)
  }
  m_pseudo = nrow(Y_pseudo)
  n_pseudo = ncol(Y_pseudo)
  
  
  M = replicate.matrix(M)
  ctl = tological(ctl, n)
  if (inputcheck) {
    # if (m > n) 
    #   warning("m is greater than n!  This is not a problem itself,
    #           but may indicate that you need to transpose your data matrix.  
    #           Please ensure that rows correspond to observations (e.g. microarrays) 
    #           and columns correspond to features (e.g. genes).")
    if (sum(is.na(Y)) > 0) 
      warning("Y contains missing values.  This is not supported.")
    if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 0) 
      warning("Y contains infinities.  This is not supported.")
  }
  
  Y_pseudo = RUV1(Y_pseudo, eta, ctl, include.intercept = include.intercept)
  
  
  
  
  if (class(BSPARAM) != "ExactParam") {
    svd_k <- k
    svd_k <- min(m_pseudo - ncol(M), sum(ctl), svd_k, na.rm = TRUE)
  } else {
    svd_k <- min(m_pseudo - ncol(M), sum(ctl), na.rm = TRUE)
  }
  
  ## m represent the number of samples/observations ncol(M)
  ## represent the number of replicates If the replicate matrix
  ## is such that we have more replicates than samples, then
  ## RUV3 is not appropriate, thus, we return the Original input
  ## matrix
  if (ncol(M) >= m_pseudo | k == 0) {
    newY <- Y_pseudo
    fullalpha <- NULL
    return(newY)
  } 
  
  if (is.null(fullalpha)) 
  {
    ## The main RUVIII process Applies the residual operator of a
    ## matrix M to a matrix Y Y0 has the same dimensions as Y,
    ## i.e. m rows (observations) and n columns (genes).
    
    Y0 <- my_residop(Y_pseudo, M)
    svdObj <- BiocSingular::runSVD(x = Y0, k = svd_k, 
                                   BPPARAM = BPPARAM, BSPARAM = BSPARAM)
    
    fullalpha <- DelayedArray::t(svdObj$u[, seq_len(svd_k), drop = FALSE]) %*% Y_pseudo
  }  ## End is.null(fullalpha)
  ###############
  alpha <- fullalpha[seq_len(min(k, nrow(fullalpha))), , drop = FALSE]
  
  alpha <- DelayedArray(alpha)
  # newY <-  Y - W %*% alpha
  
  if (normalised) {
    
    ac <- alpha[, ctl, drop = FALSE]
    
    Y <- DelayedArray(Y)
    Y_stand <- sweep(Y, 2, colMeans(Y), "-")
    
    W <- Y_stand[, ctl] %*% DelayedArray::t(ac) %*% solve(ac %*% DelayedArray::t(ac))
    W <- DelayedArray(W)
    
    
    if (!is.null(subset)) {
      newY <- Y[, subset] - W %*% alpha[, subset]
    } else {
      newY <- Y - W %*% alpha
    }
    
    if (!return.info) {
      return(newY)
    } else {
      return(list(newY = newY, M = M, fullalpha = fullalpha))
    }
  } else {
    return(list(M = M, fullalpha = fullalpha))
  }
  
}






pseudoRUVIII2 <- function(Y, Y_pseudo = NULL, Pm = NULL, ctl, k = NULL, eta = NULL, 
                          include.intercept = TRUE, inputcheck = TRUE, 
                          fullalpha = NULL, return.info = FALSE, 
                          subset = NULL,
                          BPPARAM = SerialParam(), 
                          BSPARAM = ExactParam(),
                          svd_k = 50,
                          normalised = TRUE) {
  
  if (!is.null(Y_pseudo)) {
    
    use_pseudo <- TRUE
  } else {
    use_pseudo <- FALSE
  }
  
  
  
  m = nrow(Y)
  n = ncol(Y)
  
  if (use_pseudo) {
    if (is.data.frame(Y_pseudo)) {
      Y = data.matrix(Y_pseudo)
    }
    m_pseudo = nrow(Y_pseudo)
    n_pseudo = ncol(Y_pseudo)
    
  }
  
  #M = replicate.matrix(M)
  ctl = tological(ctl, n)
  if (inputcheck) {
    if (sum(is.na(Y)) > 0) 
      warning("Y contains missing values.  This is not supported.")
    if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 
        0) 
      warning("Y contains infinities.  This is not supported.")
  }
  
  if (use_pseudo) {
    Y_pseudo = RUV1(Y_pseudo, eta, ctl, include.intercept = include.intercept)
  } else {
    Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
  }
  
  
  
  if (use_pseudo) {
    if (class(BSPARAM) != "ExactParam") {
      #svd_k <- min(m_pseudo - ncol(M), sum(ctl), svd_k, na.rm = TRUE)
      svd_k <- min(m_pseudo, sum(ctl), svd_k, na.rm = TRUE)
    } else {
      #svd_k <- min(m_pseudo - ncol(M), sum(ctl), na.rm = TRUE)
      svd_k <- min(m_pseudo, sum(ctl), na.rm = TRUE)
    }
  } else {
    if (class(BSPARAM) != "ExactParam") {
      #svd_k <- min(m_pseudo - ncol(M), sum(ctl), svd_k, na.rm = TRUE)
      svd_k <- min(m, sum(ctl), svd_k, na.rm = TRUE)
    } else {
      #svd_k <- min(m_pseudo - ncol(M), sum(ctl), na.rm = TRUE)
      svd_k <- min(m, sum(ctl), na.rm = TRUE)
    }
  }
  
  
  
  ## m represent the number of samples/observations ncol(M)
  ## represent the number of replicates If the replicate matrix
  ## is such that we have more replicates than samples, then
  ## RUV3 is not appropriate, thus, we return the Original input
  ## matrix
  # if (ncol(M) >= m_pseudo | k == 0) {
  #   newY <- Y_pseudo
  #   fullalpha <- NULL
  #   return(newY)
  # } 
  
  if (is.null(fullalpha)) {
    ## The main RUVIII process Applies the residual operator of a
    ## matrix M to a matrix Y Y0 has the same dimensions as Y,
    ## i.e. m rows (observations) and n columns (genes).
    
    #Y0 <- my_residop(Y_pseudo, M)
    
    if (use_pseudo) {
      Y0 <- residop2(Pm, Y_pseudo)
      svdObj <- BiocSingular::runSVD(x = Y0, k = svd_k,
                                     BPPARAM = BPPARAM, 
                                     BSPARAM = BSPARAM)
      
      fullalpha <- DelayedArray::t(svdObj$u[, seq_len(svd_k), drop = FALSE]) %*% Y_pseudo
    } else {
      Y0 <- residop2(Pm, Y)
      svdObj <- BiocSingular::runSVD(x = Y0, 
                                     k = svd_k, 
                                     BPPARAM = BPPARAM, 
                                     BSPARAM = BSPARAM)
      
      fullalpha <- DelayedArray::t(svdObj$u[, seq_len(svd_k), drop = FALSE]) %*% Y
    }
    
  }  ## End is.null(fullalpha)
  ###############
  alpha <- fullalpha[seq_len(min(k, nrow(fullalpha))), , drop = FALSE]
  ac <- alpha[, ctl, drop = FALSE]
  
  Y <- DelayedArray(Y)
  Y_stand <- sweep(Y, 2, colMeans(Y), "-")
  
  W <- Y_stand[, ctl] %*% DelayedArray::t(ac) %*% solve(ac %*% DelayedArray::t(ac))
  W <- DelayedArray(W)
  alpha <- DelayedArray(alpha)
  # newY <-  Y - W %*% alpha
  
  if (normalised) {
    if (!is.null(subset)) {
      newY <- Y[, subset] - W %*% alpha[, subset]
    } else {
      newY <- Y - W %*% alpha
    }
    
    if (!return.info) {
      return(newY)
    } else {
      return(list(newY = newY, Pm = Pm, fullalpha = fullalpha))
    }
  } else {
    return(list(Pm = Pm, fullalpha = fullalpha))
  }
  
}




# 
# 
# pseudoRUVIII3 <-  function(Y, Y_pseudo, M, ctl, k = NULL, eta = NULL, include.intercept = TRUE, batch = NULL,
#                            average = FALSE, fullalpha = NULL, return.info = FALSE, inputcheck = TRUE,
#                            regularise = 0.01, scale = TRUE, 
#                            BPPARAM = SerialParam(), BSPARAM = ExactParam()) 
# {
#   # if (scale) {
#   #   scale_res <- standardize2(t(Y), batch)
#   #   stand_tY <- DelayedArray::t(scale_res$stand_Y)
#   #   stand_sd <- sqrt(scale_res$stand_var)
#   #   stand_mean <- scale_res$stand_mean
#   #   Y <- stand_tY
#   #   
#   #   rm("stand_tY")
#   # } 
#   
#   if (is.data.frame(Y)) 
#     Y = data.matrix(Y)
#   m = nrow(Y)
#   n = ncol(Y)
#   
#   if (is.data.frame(Y_pseudo)) 
#     Y = data.matrix(Y_pseudo)
#   m_pseudo = nrow(Y_pseudo)
#   n_pseudo = ncol(Y_pseudo)
#   
#   
#   M = replicate.matrix(M)
#   ctl = tological(ctl, n)
#   if (inputcheck) {
#     if (m > n) 
#       warning("m is greater than n!  This is not a problem itself, but may indicate that you need to transpose your data matrix.  Please ensure that rows correspond to observations (e.g. microarrays) and columns correspond to features (e.g. genes).")
#     if (sum(is.na(Y)) > 0) 
#       warning("Y contains missing values.  This is not supported.")
#     if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 
#         0) 
#       warning("Y contains infinities.  This is not supported.")
#   }
#   Y_pseudo = RUV1(Y_pseudo, eta, ctl, include.intercept = include.intercept)
#   
#   
#   if (class(BSPARAM) != "ExactParam") {
#     svd_k <- min(m_pseudo - ncol(M), sum(ctl), svd_k, na.rm = TRUE)
#   } else {
#     svd_k <- min(m_pseudo - ncol(M), sum(ctl), na.rm = TRUE)
#   }
#   
#   ## m represent the number of samples/observations ncol(M)
#   ## represent the number of replicates If the replicate matrix
#   ## is such that we have more replicates than samples, then
#   ## RUV3 is not appropriate, thus, we return the Original input
#   ## matrix
#   if (ncol(M) >= m_pseudo | k == 0) {
#     newY <- Y_pseudo
#     fullalpha <- NULL
#   } else {
#     
#     if (is.null(fullalpha)) 
#     {
#       ## The main RUVIII process Applies the residual operator of a
#       ## matrix M to a matrix Y Y0 has the same dimensions as Y,
#       ## i.e. m rows (observations) and n columns (genes).
#       
#       # if(class(Y) == "matrix"){
#       #   Y0 <- eigenResidop(Y_pseudo, M)
#       # } else if (class(Y) == "dgeMatrix"){
#       #   Y0 <- eigenResidop(as.matrix(Y_pseudo), M)
#       # } else {
#       #   Y0 <- my_residop(Y_pseudo, M)
#       # }
#       
#       Y0 <- my_residop(Y_pseudo, M)
#       svdObj <- BiocSingular::runSVD(x = Y0, k = svd_k, BPPARAM = BPPARAM, BSPARAM = BSPARAM)
#       
#       fullalpha <- t(svdObj$u[, seq_len(svd_k), drop = FALSE]) %*% Y_pseudo
#     }  ## End is.null(fullalpha)
#     ###############
#     alpha <- fullalpha[seq_len(min(k, nrow(fullalpha))), , drop = FALSE]
#     ac <- alpha[, ctl, drop = FALSE]
#     W <- Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
#     print(dim(alpha))
#     print(dim(W))
#     
#     
#     # newY <- Y - W %*% alpha
#     
#   }  ## End else(ncol(M) >= m | k == 0)
#   
#   rm("Y_pseudo")
#   
#   
#   return(list(newY = Y, W = W, fullalpha = alpha))
# }
# 


########################################################################################


findMNC_old <- function(exprs_mat, clustering_list, dist = "euclidean", 
                        BPPARAM, plot_igraph = TRUE) {
  
  batch_num <- length(clustering_list)
  names(clustering_list) <- paste("Batch", seq_len(batch_num), 
                                  sep = "")
  
  ## Check which batch has only one cluster
  batch_oneType <- which(unlist(lapply(clustering_list, function(x) length(levels(as.factor(x))) == 1)))
  ## If there existsome batch_oneType
  if (length(batch_oneType) != 0) {
    ## And if all batch_oneType == num of batches, i.e.  every
    ## batch only contains one cell type
    if (length(batch_oneType) == batch_num) {
      combine_pair <- utils::combn(batch_num, 2)
      batch_oneType <- NULL
      allones <- TRUE
    } else {
      ## if at least some batch contains more than 1 cell type Then
      ## take away the batches with only one cell type and then
      ## iterate through all combn
      ## 
      
      ## 20190606 YL: to prevent only one batch left after exclude the batch with one type.
      ## 
      if (length(c(seq_len(batch_num))[-batch_oneType]) == 1) {
        combine_pair <- NULL
      }else{
        combine_pair <- utils::combn(c(seq_len(batch_num))[-batch_oneType], 
                                     2)
      }
      
      
      ## And then for those batches with only one cell type, we bind
      ## to the previous generated combn
      for (i in batch_oneType) {
        for (j in c(seq_len(batch_num))[-batch_oneType]) {
          combine_pair <- cbind(combine_pair, c(i, j))
        }
      }
      allones <- FALSE
    }
  } else {
    combine_pair <- utils::combn(batch_num, 2)
    allones <- FALSE
  }
  
  
  
  
  mnc <- list()
  
  ## If there are only two batches containing only two cell
  ## types, then finding MNN is trivial. Return NULL
  if (allones & batch_num == 2) {
    return(NULL)
  } else if (allones) {
    ## If every batch contains only one cell type...
    dist_res <- matrix(NA, nrow = batch_num, ncol = batch_num)
    
    
    dist_mat_med <- BiocParallel::bplapply(seq_len(ncol(combine_pair)), function(k) {
      compute_dist_mat_med(k = k, exprs_mat = exprs_mat, 
                           clustering_list = clustering_list, 
                           combine_pair = combine_pair, 
                           dist = dist)
    }, BPPARAM = BPPARAM)
    
    
    
    for (k in seq_len(ncol(combine_pair))) {
      dist_res[combine_pair[1, k], combine_pair[2, k]] <- dist_res[combine_pair[2, k], combine_pair[1, k]] <- dist_mat_med[[k]]
    }
    
    ## The neighbour_res is then which ever two pairs of batches
    ## that are closes to each other
    neighbour_res <- apply(dist_res, 1, which.min)
    
    mnc_mat <- c()
    for (i in seq_along(neighbour_res)) {
      if (neighbour_res[neighbour_res[i]] == i) {
        mnc_mat <- rbind(mnc_mat, sort(c(i, neighbour_res[i])))
      }
    }
    mnc_mat <- unique(mnc_mat)
    mnc <- list()
    for (i in seq_len(nrow(mnc_mat))) {
      mnc[[i]] <- matrix(1, ncol = 2, nrow = 1)
      colnames(mnc[[i]]) <- c(paste("Batch", mnc_mat[i, 
                                                     1], sep = ""), paste("Batch", mnc_mat[i, 2], 
                                                                          sep = ""))
    }
  } else {
    for (k in seq_len(ncol(combine_pair))) {
      dist_res <- list()
      # print(k)
      res1 <- clustering_list[[combine_pair[1, k]]]
      res2 <- clustering_list[[combine_pair[2, k]]]
      
      
      dist_res <- BiocParallel::bplapply(seq_len(max(res1)), function(i) {
        compute_dist_res(i = i, res1 = res1, res2 = res2,
                         exprs_mat = exprs_mat, dist = dist, dist_res)
      }, BPPARAM = BPPARAM)
      
      
      dist_res <- do.call(rbind, dist_res)
      neighbour_batch1 <- apply(dist_res, 1, which.min)
      neighbour_batch2 <- apply(dist_res, 2, which.min)
      mnc_tmp <- c()
      for (l in seq_len(length(neighbour_batch1))) {
        if (neighbour_batch2[neighbour_batch1[l]] == l) {
          mnc_tmp <- rbind(mnc_tmp, c(l, neighbour_batch1[l]))
        }
      }
      mnc[[k]] <- mnc_tmp
      colnames(mnc[[k]]) <- c(paste("Batch", combine_pair[1, k], sep = ""), 
                              paste("Batch", combine_pair[2, k], sep = ""))
    }
  }  ## End else
  
  ############################################################### Function to perform network analysis
  
  
  edge_list <- do.call(rbind, lapply(mnc, function(x) t(apply(x, 1, function(y) paste(colnames(x), y, sep = "_")))))
  
  
  ### 20190606 YL: creating the node list for each cluster in each batch
  
  node_list <- list()
  
  for (b in seq_along(clustering_list)) {
    node_list[[b]] <- paste(names(clustering_list)[b],
                            seq_len(max(clustering_list[[b]])), sep = "_")
  }
  
  node_list <- unlist(node_list)
  
  
  
  ### 20190606 YL: Use the node list and edge list to build the network
  ### The edge list is null case will not be a special case now
  
  
  g <- igraph::make_empty_graph(directed = FALSE) + igraph::vertices(node_list)
  g <- g + igraph::edges(t(edge_list))
  
  
  if(plot_igraph){
    igraph::plot.igraph(g)
  }
  
  mnc <- igraph::fastgreedy.community(g)
  mnc_df <- data.frame(group = as.numeric(mnc$membership),
                       batch = as.numeric(gsub("Batch", "", gsub("_.*", "", mnc$names))),
                       cluster = as.numeric(gsub(".*_", "", mnc$names)))
  
  return(mnc_df)
}

compute_dist_mat_med <- function(k, exprs_mat, clustering_list, 
                                 combine_pair, dist) {
  ## We go through every pairwise batches print(k) Extract the
  ## cell type information and the expression matrices
  res1 <- clustering_list[[combine_pair[1, k]]]
  res2 <- clustering_list[[combine_pair[2, k]]]
  
  mat1 <- as(exprs_mat[, names(which(res1 == 1))], "dgCMatrix")
  mat2 <- as(exprs_mat[, names(which(res2 == 1))], "dgCMatrix")
  
  if (dist == "cosine") {
    dist_mat <- 1 - proxyC::simil(t(mat1), t(mat2), method = "cosine")
  } else if (dist == "cor") {
    dist_mat <- 1 - proxyC::simil(t(mat1), t(mat2), method = "correlation")
  } else {
    dist_mat <- proxyC::dist(t(mat1), t(mat2))
  }
  
  
  dist_mat <- as.matrix(dist_mat)
  
  ## The dist_res (distance measure between batches) is then the
  ## median of all pairwise distances
  return(stats::median(dist_mat))
}

findMNC <- function(exprs_mat, clustering_list, dist = "euclidean", 
                    BPPARAM, plot_igraph = TRUE) {
  
  batch_num <- length(clustering_list)
  names(clustering_list) <- paste("Batch", seq_len(batch_num), sep = "")
  
  ## Check which batch has only one cluster
  batch_oneType <- which(unlist(lapply(clustering_list, function(x) length(levels(as.factor(x))) == 1)))
  ## If there existsome batch_oneType
  if (length(batch_oneType) != 0) {
    ## And if all batch_oneType == num of batches, i.e.  every
    ## batch only contains one cell type
    if (length(batch_oneType) == batch_num) {
      combine_pair <- utils::combn(batch_num, 2)
      batch_oneType <- NULL
      allones <- TRUE
    } else {
      ## if at least some batch contains more than 1 cell type Then
      ## take away the batches with only one cell type and then
      ## iterate through all combn
      ## 
      
      ## 20190606 YL: to prevent only one batch left after exclude the batch with one type.
      ## 
      if (length(c(seq_len(batch_num))[-batch_oneType]) == 1) {
        combine_pair <- NULL
      }else{
        combine_pair <- utils::combn(c(seq_len(batch_num))[-batch_oneType], 2)
      }
      
      
      ## And then for those batches with only one cell type, we bind
      ## to the previous generated combn
      for (i in batch_oneType) {
        for (j in c(seq_len(batch_num))[-batch_oneType]) {
          combine_pair <- cbind(combine_pair, c(i, j))
        }
      }
      allones <- FALSE
    }
  } else {
    combine_pair <- utils::combn(batch_num, 2)
    allones <- FALSE
  }
  
  
  
  
  mnc <- list()
  
  ## If there are only two batches containing only two cell
  ## types, then finding MNN is trivial. Return NULL
  if (allones & batch_num == 2) {
    return(NULL)
  } else if (allones) {
    ## If every batch contains only one cell type...
    dist_res <- matrix(NA, nrow = batch_num, ncol = batch_num)
    
    
    dist_mat_med <- BiocParallel::bplapply(seq_len(ncol(combine_pair)), function(k) {
      compute_dist_mat_med(k = k, exprs_mat = exprs_mat, 
                           clustering_list = clustering_list, 
                           combine_pair = combine_pair, 
                           dist = dist)
    }, BPPARAM = BPPARAM)
    
    
    
    for (k in seq_len(ncol(combine_pair))) {
      dist_res[combine_pair[1, k], combine_pair[2, k]] <- 
        dist_res[combine_pair[2, k], combine_pair[1, k]] <- 
        dist_mat_med[[k]]
    }
    
    ## The neighbour_res is then which ever two pairs of batches
    ## that are closes to each other
    neighbour_res <- apply(dist_res, 1, which.min)
    
    mnc_mat <- c()
    for (i in seq_along(neighbour_res)) {
      if (neighbour_res[neighbour_res[i]] == i) {
        mnc_mat <- rbind(mnc_mat, sort(c(i, neighbour_res[i])))
      }
    }
    mnc_mat <- unique(mnc_mat)
    mnc <- list()
    for (i in seq_len(nrow(mnc_mat))) {
      mnc[[i]] <- matrix(1, ncol = 2, nrow = 1)
      colnames(mnc[[i]]) <- c(paste("Batch", mnc_mat[i, 1], sep = ""),
                              paste("Batch", mnc_mat[i, 2], sep = ""))
    }
  } else {
    
    
    mnc <- pbmcapply::pbmclapply(seq_len(ncol(combine_pair)), function(k) {
      res1 <- clustering_list[[combine_pair[1, k]]]
      res2 <- clustering_list[[combine_pair[2, k]]]
      
      mat1 <- as(exprs_mat[, names(res1), drop = FALSE], "dgCMatrix")
      mat2 <- as(exprs_mat[, names(res2), drop = FALSE], "dgCMatrix")
      
      if (dist == "cosine") {
        dist_mat <- 1 - proxyC::simil(t(mat1), t(mat2), method = "cosine")
      } else if (dist == "cor") {
        dist_mat <- 1 - proxyC::simil(t(mat1), t(mat2), method = "correlation")
      } else {
        dist_mat <- proxyC::dist(t(mat1), t(mat2))
      }
      
      dist_res <- apply(expand.grid(seq_len(max(res1)), seq_len(max(res2))), 1, 
                        function(x) {
                          median(dist_mat[res1 == x[1], res2 == x[2], drop = FALSE], na.rm = TRUE)
                        })
      
      dist_res <- matrix(dist_res, nrow = (max(res1)), ncol = (max(res2)))
      
      neighbour_batch1 <- apply(dist_res, 1, which.min)
      neighbour_batch2 <- apply(dist_res, 2, which.min)
      mnc_tmp <- c()
      for (l in seq_len(length(neighbour_batch1))) {
        if (neighbour_batch2[neighbour_batch1[l]] == l) {
          mnc_tmp <- rbind(mnc_tmp, c(l, neighbour_batch1[l]))
        }
      }
      colnames(mnc_tmp) <- c(paste("Batch", combine_pair[1, k], sep = ""), 
                             paste("Batch", combine_pair[2, k], sep = ""))
      return(mnc_tmp)
    }, mc.cores = 5)
    
    
  }  ## End else
  
  ############################################################### Function to perform network analysis
  
  
  edge_list <- do.call(rbind, lapply(mnc, function(x) t(apply(x, 1, function(y) paste(colnames(x), y, sep = "_")))))
  
  
  ### 20190606 YL: creating the node list for each cluster in each batch
  
  node_list <- list()
  
  for (b in seq_along(clustering_list)) {
    node_list[[b]] <- paste(names(clustering_list)[b],
                            seq_len(max(clustering_list[[b]])), sep = "_")
  }
  
  node_list <- unlist(node_list)
  
  
  
  ### 20190606 YL: Use the node list and edge list to build the network
  ### The edge list is null case will not be a special case now
  
  
  g <- igraph::make_empty_graph(directed = FALSE) + igraph::vertices(node_list)
  g <- g + igraph::edges(t(edge_list))
  
  
  if(plot_igraph){
    igraph::plot.igraph(g)
  }
  
  mnc <- igraph::fastgreedy.community(g)
  mnc_df <- data.frame(group = as.numeric(mnc$membership),
                       batch = as.numeric(gsub("Batch", "", gsub("_.*", "", mnc$names))),
                       cluster = as.numeric(gsub(".*_", "", mnc$names)))
  
  return(mnc_df)
}



compute_dist_mat_med <- function(k, exprs_mat, clustering_list, 
                                 combine_pair, dist) {
  ## We go through every pairwise batches print(k) Extract the
  ## cell type information and the expression matrices
  res1 <- clustering_list[[combine_pair[1, k]]]
  res2 <- clustering_list[[combine_pair[2, k]]]
  mat1 <- exprs_mat[, names(which(res1 == 1))]
  mat2 <- exprs_mat[, names(which(res2 == 1))]
  
  ## The distance between matrices are calculated as such...
  if (dist == "cosine") {
    dist_mat <- proxy::dist(t(mat1), t(mat2), method = "cosine")
  } else if (dist == "cor") {
    dist_mat <- 1 - stats::cor(as.matrix(mat1), as.matrix(mat2))
  } else {
    dist_mat <- pdist::pdist(t(mat1), t(mat2))
  }
  
  
  dist_mat <- as.matrix(dist_mat)
  
  ## The dist_res (distance measure between batches) is then the
  ## median of all pairwise distances
  return(stats::median(dist_mat))
}



###################################################################################################### Function to create replicate based on the mutual nearest
###################################################################################################### cluster results
#### Note this function mnc_df is the output from findMNC function

mncReplicate <- function(clustering_list, clustering_distProp, 
                         replicate_prop, mnc_df) {
  # batch_num <- length(clustering_list) ### 20190606 YL... this line is not used in this function
  
  ### 20190606 YL: NULL is still necessary for case allones & two batches 
  if (!is.null(mnc_df)) {
    
    replicate_vector <- rep(NA, length(unlist(clustering_list)))
    names(replicate_vector) <- names(unlist(clustering_list))
    clustering_distProp <- unlist(clustering_distProp)
    replicate_size <- table(mnc_df$group)
    for (i in seq_len(max(mnc_df$group))) {
      tmp_names <- c()
      
      # For each batch in this replicate
      mnc_df_sub <- mnc_df[mnc_df$group == i, ]
      for (l in seq_len(replicate_size[i])) {
        
        tmp_names <- c(tmp_names, 
                       names(which(clustering_list[[mnc_df_sub[l, "batch"]]] == 
                                     mnc_df_sub[l, "cluster"])))
      }
      
      replicate_vector[tmp_names[clustering_distProp[tmp_names] <= 
                                   replicate_prop]] <- paste("Replicate", i, sep = "_")
    }
    
    
    replicate_vector[is.na(replicate_vector)] <- seq_len(sum(is.na(replicate_vector)))
  } else {
    replicate_vector <- rep(NA, length(unlist(clustering_list)))
    names(replicate_vector) <- names(unlist(clustering_list))
    clustering_distProp <- unlist(clustering_distProp)
    current_idx <- 1
    for (j in seq_along(clustering_list)) {
      for (k in seq_len(max(clustering_list[[j]]))) {
        tmp_names <- names(which(clustering_list[[j]] == 
                                   k))
        replicate_vector[tmp_names[clustering_distProp[tmp_names] <= 
                                     replicate_prop]] <- paste("Replicate", current_idx + 
                                                                 k, sep = "_")
      }
      current_idx <- current_idx + max(clustering_list[[j]])
    }
    replicate_vector[is.na(replicate_vector)] <- seq_len(sum(is.na(replicate_vector)))
  }
  
  return(replicate_vector)
}

# 
# 
# compute_dist_res2 <- function(i, res1, res2, exprs_mat, dist, 
#                               dist_res) {
#   res_tmp <- c()
#   mat1 <- as(exprs_mat[, names(which(res1 == i)), drop = FALSE], "dgCMatrix")
#   for (j in seq_len(max(res2))) {
#     mat2 <- as(exprs_mat[, names(which(res2 == j)), drop = FALSE], "dgCMatrix")
#     if (dist == "cosine") {
#       #dist_mat <- stats::dist(t(mat1), t(mat2), method = "cosine")
#       dist_mat <- 1 - proxyC::simil(t(mat1), t(mat2), method = "cosine")
#     } else if (dist == "cor") {
#       #dist_mat <- 1 - stats::cor(as.matrix(mat1), as.matrix(mat2))
#       dist_mat <- 1 - proxyC::simil(t(mat1), t(mat2), method = "correlation")
#     } else {
#       #dist_mat <- pdist::pdist(t(mat1), t(mat2))
#       dist_mat <- proxyC::dist(t(mat1), t(mat2))
#     }
#     dist_mat <- as.matrix(dist_mat)
#     res_tmp <- c(res_tmp, stats::median(dist_mat))
#   }
#   dist_res[[i]] <- res_tmp
#   
#   return(dist_res[[i]])
# }



compute_dist_res2 <- function(i, res1, res2, exprs_mat, dist, 
                              dist_res) {
  res_tmp <- c()
  mat1 <- as(exprs_mat[, names(which(res1 == i)), drop = FALSE], "dgCMatrix")
  mat2 <- as(exprs_mat[, names(res2), drop = FALSE], "dgCMatrix")
  
  if (dist == "cosine") {
    #dist_mat <- stats::dist(t(mat1), t(mat2), method = "cosine")
    dist_mat <- 1 - proxyC::simil(t(mat1), t(mat2), method = "cosine")
  } else if (dist == "cor") {
    #dist_mat <- 1 - stats::cor(as.matrix(mat1), as.matrix(mat2))
    dist_mat <- 1 - proxyC::simil(t(mat1), t(mat2), method = "correlation")
  } else {
    #dist_mat <- pdist::pdist(t(mat1), t(mat2))
    dist_mat <- proxyC::dist(t(mat1), t(mat2))
  }
  
  dist_res[[i]] <- sapply(seq_len(max(res2)), function(x) median(dist_mat[, res2 == x]))
  
  return(dist_res[[i]])
}




compute_dist_res <- function(i, res1, res2, exprs_mat, dist, 
                             dist_res) {
  res_tmp <- c()
  mat1 <- exprs_mat[, names(which(res1 == i))]
  for (j in seq_len(max(res2))) {
    mat2 <- exprs_mat[, names(which(res2 == j))]
    if (dist == "cosine") {
      dist_mat <- stats::dist(t(mat1), t(mat2), method = "cosine")
    } else if (dist == "cor") {
      dist_mat <- 1 - stats::cor(as.matrix(mat1), as.matrix(mat2))
    } else {
      dist_mat <- pdist::pdist(t(mat1), t(mat2))
    }
    dist_mat <- as.matrix(dist_mat)
    res_tmp <- c(res_tmp, stats::median(dist_mat))
  }
  dist_res[[i]] <- res_tmp
  
  return(dist_res[[i]])
}



estY <- function(Y, fullalpha, ctl = colnames(Y), subset = NULL, k = 20) {
  
  if (class(ctl) %in% c("character")) {
    ctl <- intersect(colnames(Y), ctl)
    if (length(ctl) == 0) {
      stop("No provided ctl genes in the data")
    }
  }
  
  m <- nrow(Y)
  Y <- DelayedArray(Y)
  Y_stand <- sweep(Y, 2, DelayedMatrixStats::colMeans2(Y), "-")
  fullalpha <- DelayedArray(fullalpha)
  alpha <- fullalpha[seq_len(min(k, nrow(fullalpha))), , drop = FALSE]
  W <- Y_stand[, ctl] %*% DelayedArray::t(alpha[, ctl]) %*% solve(alpha[, ctl] %*% DelayedArray::t(alpha[, ctl]))
  
  if (!is.null(subset)) {
    newY <- Y[, subset] - W %*% alpha[, subset]
  } else {
    newY <- Y - W %*% alpha
  }
  
  return(newY)
}


estY2 <- function(Y, fullalpha, ctl = colnames(Y), subset = NULL) {
  
  m <- nrow(Y)
  # Y <- DelayedArray(Y)
  Y_stand <- sweep(Y, 2, colMeans(Y), "-")
  #  fullalpha <- DelayedArray(fullalpha)
  ac <- fullalpha[, ctl]
  W <- Y_stand[, ctl] %*% t(fullalpha[, ctl]) %*% solve(fullalpha[, ctl] %*% t(fullalpha[, ctl]))
  
  if (!is.null(subset)) {
    newY <- Y[, subset] - W %*% fullalpha[, subset]
  } else {
    newY <- Y - W %*% fullalpha
  }
  
  return(newY)
}



library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(scattermore)
library(SingleCellExperiment)
library(scran)
library(BiocParallel)


#source("scMerge2_functions/feature_selection.R")

library(scMerge)
data("segList")
library(ruv)
library(BiocParallel)
library(DelayedArray)

# Settings

feature_selection <- function(exprsMat, batch, top_n = 2000, BPPARAM = SerialParam()) {
  combined.dec <- scran::modelGeneVar(exprsMat, block = batch, BPPARAM = BPPARAM)
  chosen.hvgs <- scran::getTopHVGs(combined.dec, n = top_n)
  return(chosen.hvgs)
}



construct_pseudoBulk_sample <- function(i,
                                        exprsMat,
                                        clust,
                                        pseudoBulk_fn,
                                        pseudobulk_sample,
                                        pseudobulk_sample_list,
                                        exprsMat_counts,
                                        pseudoBulk_arg,
                                        k_psuedoBulk,
                                        ncores) {
  
  min_cells <- min(min(table(pseudobulk_sample)), pseudoBulk_arg$min_cells)
  
  
  if (as.character(substitute(pseudoBulk_fn)) %in% c("create_pseudoBulk_pool", 
                                                     "create_pseudoBulk_pool_divide",
                                                     "create_pseudoBulk_divide")) {
    res<- pseudoBulk_fn(exprsMat_counts[, pseudobulk_sample == pseudobulk_sample_list[i]],
                        clust[[i]], k_fold = k_psuedoBulk,
                        ncore = ncores)
  } else if (as.character(substitute(pseudoBulk_fn)) %in% "create_pseudoBulk_graph_divide") {
    res <- pseudoBulk_fn(exprsMat[, pseudobulk_sample == pseudobulk_sample_list[i]],
                         exprsMat_counts[, pseudobulk_sample == pseudobulk_sample_list[i]],
                         clust[[i]], k_fold = k_psuedoBulk,
                         ncore = ncores)
    
  } else if (as.character(substitute(pseudoBulk_fn)) %in% "create_pseudoBulk_KNN") {
    
    res <- pseudoBulk_fn(exprsMat[, pseudobulk_sample == pseudobulk_sample_list[i]], 
                         clust[[i]], 
                         prop_pseudoBulk = pseudoBulk_arg$prop_pseudoBulk,
                         chosen.hvg = pseudoBulk_arg$chosen.hvg, 
                         pseudobulk_k = pseudoBulk_arg$pseudobulk_k,
                         upsample = pseudoBulk_arg$upsample,
                         min_cells = min_cells) 
    
  } else if (as.character(substitute(pseudoBulk_fn)) %in% "create_pseudoBulk_sample") {
    
    res <- pseudoBulk_fn(exprsMat[, pseudobulk_sample == pseudobulk_sample_list[i]], 
                         clust[[i]], 
                         prop_pseudoBulk = pseudoBulk_arg$prop_pseudoBulk,
                         ncore = ncores,
                         upsample = pseudoBulk_arg$upsample,
                         min_cells  = min_cells) 
    
  } else {
    res <- pseudoBulk_fn(exprsMat[, pseudobulk_sample == pseudobulk_sample_list[i]],
                         clust[[i]], k_fold = k_psuedoBulk,
                         ncore = ncores)
    
  }
  rownames(res) <- paste(i, rownames(res), sep = "_")
  
  
  
  return(res)
  
}

scMerge2 <- function(exprsMat, 
                     batch, 
                     cellTypes = NULL, 
                     condition = NULL,
                     use_bpparam = BiocParallel::SerialParam(),
                     use_bsparam = BiocSingular::ExactParam(),
                     use_bnparam = BiocNeighbors::AnnoyParam(),
                     ruvK = 30, 
                     ctl = rownames(exprsMat), 
                     sample = NULL,
                     cluster_by_sample = FALSE,
                     match_within_batch = TRUE,
                     match_by_sample = TRUE,    
                     exprsMat_counts = NULL,
                     k_psuedoBulk = 30, 
                     k_celltype = 5,
                     pseudoBulk_fn = create_pseudoBulk_graph,
                     pseudoBulk_arg = list(),
                     ncores = 1,
                     chosen.hvg = NULL,
                     cosineNorm = TRUE,
                     return_subset = FALSE,
                     normalised = TRUE) {
  
  batch <- as.character(batch)
  
  
  batch_list <- sort(unique(batch))
  
  
  if (!is.null(sample)) {
    sample_mode <- TRUE
    condition_mode <- FALSE
    sample <- as.character(sample)
    # meta_sample <- cbind(batch = batch, sample = sample)
    # meta_sample <- unique(meta_sample)
    sample_list <- sort(unique(sample))
    pseudobulk_sample <- sample
    pseudobulk_sample_list <- sample_list
    pseudobulk_batch <- batch
    pseudobulk_batch_list <- batch_list
  } else {
    sample_mode <- FALSE
    pseudobulk_sample <- batch
    pseudobulk_sample_list <- batch_list
    pseudobulk_batch <- batch
    pseudobulk_batch_list <- batch_list
  }
  
  
  
  if (!is.null(condition)) {
    cond_mode <- TRUE
    sample_mode <- FALSE
    condition_list <- sort(unique(condition))
  } else {
    cond_mode <- FALSE
  }
  
  meta_sample <- cbind(batch = batch, sample = sample, condition = condition)
  meta_sample <- unique(meta_sample)
  
  
  
  if (is.null(chosen.hvg)) {
    chosen.hvg <- feature_selection(exprsMat, batch, BPPARAM = use_bpparam)
  }
  
  
  if (is.null(cellTypes)) {
    
    if (cluster_by_sample) {
      print("Cluster by Sample")
      clust <- BiocParallel::bplapply(pseudobulk_sample_list, function(x) {
        g <- scran::buildSNNGraph(exprsMat[chosen.hvg, pseudobulk_sample == x], 
                                  k = k_celltype, 
                                  BNPARAM = use_bnparam)
        igraph::cluster_louvain(g)$membership
        
      }, BPPARAM = use_bpparam)
    } else {
      print("Cluster by Batch")
      clust <- BiocParallel::bplapply(pseudobulk_batch_list, function(x) {
        g <- scran::buildSNNGraph(exprsMat[chosen.hvg, pseudobulk_batch == x], 
                                  k = k_celltype, 
                                  BNPARAM = use_bnparam)
        res <- igraph::cluster_louvain(g)$membership
        names(res) <- colnames(exprsMat)[pseudobulk_batch == x]
        res
      }, BPPARAM = use_bpparam)
      clust <- unlist(clust)
      clust <- clust[colnames(exprsMat)]
      clust <- split(clust, pseudobulk_sample)
      clust <- lapply(clust, function(x) {
        as.numeric(factor(x))
      })
      print(lapply(clust, table))
    }
    
  } else {
    clust <- BiocParallel::bplapply(pseudobulk_sample_list, function(x) {
      res <- as.numeric(factor(cellTypes[pseudobulk_sample == x]))
      names(res) <- colnames(exprsMat)[pseudobulk_sample == x]
      res
    }, BPPARAM = use_bpparam)
  }
  
  
  if (cosineNorm) {
    print("Normalising data")
    exprsMat <- batchelor::cosineNorm(exprsMat, BPPARAM = use_bpparam)
  }
  
  
  
  
  
  #set.seed(2020)
  
  print("Constructing pseudo-bulk")
  pseudoBulk_arg$chosen.hvg <- chosen.hvg
  min_cells <- min(min(table(pseudobulk_sample)), pseudoBulk_arg$min_cells)
  
  
  aggExprs <- list()
  for( i in seq_along(pseudobulk_sample_list)){
    print(i)
    
    if (as.character(substitute(pseudoBulk_fn)) %in% c("create_pseudoBulk_pool", 
                                                       "create_pseudoBulk_pool_divide",
                                                       "create_pseudoBulk_divide")) {
      res<- pseudoBulk_fn(exprsMat_counts[, pseudobulk_sample == pseudobulk_sample_list[i]],
                          clust[[i]], k_fold = k_psuedoBulk,
                          ncore = ncores)
    } else if (as.character(substitute(pseudoBulk_fn)) %in% "create_pseudoBulk_graph_divide") {
      res <- pseudoBulk_fn(exprsMat[, pseudobulk_sample == pseudobulk_sample_list[i]],
                           exprsMat_counts[, pseudobulk_sample == pseudobulk_sample_list[i]],
                           clust[[i]], k_fold = k_psuedoBulk,
                           ncore = ncores)
      
    } else if (as.character(substitute(pseudoBulk_fn)) %in% "create_pseudoBulk_KNN") {
      
      res <- pseudoBulk_fn(exprsMat[, pseudobulk_sample == pseudobulk_sample_list[i]], 
                           clust[[i]], 
                           prop_pseudoBulk = pseudoBulk_arg$prop_pseudoBulk,
                           chosen.hvg = pseudoBulk_arg$chosen.hvg, 
                           pseudobulk_k = pseudoBulk_arg$pseudobulk_k,
                           upsample = pseudoBulk_arg$upsample,
                           min_cells = min_cells) 
      
    } else if (as.character(substitute(pseudoBulk_fn)) %in% "create_pseudoBulk_sample") {
      
      res <- pseudoBulk_fn(exprsMat[, pseudobulk_sample == pseudobulk_sample_list[i]], 
                           clust[[i]], 
                           prop_pseudoBulk = pseudoBulk_arg$prop_pseudoBulk,
                           ncore = ncores,
                           upsample = pseudoBulk_arg$upsample,
                           min_cells  = min_cells) 
      
    } else {
      res <- pseudoBulk_fn(exprsMat[, pseudobulk_sample == pseudobulk_sample_list[i]],
                           clust[[i]], k_fold = k_psuedoBulk,
                           ncore = ncores)
      
    }
    rownames(res) <- paste(i, rownames(res), sep = "_")
    print(dim(res))
    aggExprs[[i]] <-res
  }
  
  
  
  
  bulkExprs <- t(do.call(rbind, aggExprs))
  if (as.character(substitute(pseudoBulk_fn)) %in% c("create_pseudoBulk_pool",
                                                     "create_pseudoBulk_pool_divide",
                                                     "create_pseudoBulk_graph_divide",
                                                     "create_pseudoBulk_divide") & cosineNorm) {
    bulkExprs <- batchelor::cosineNorm(bulkExprs, BPPARAM = use_bpparam)
  }
  
  
  rownames(bulkExprs) <- rownames(exprsMat)
  
  
  pseudo_batch <- unlist(lapply(strsplit(colnames(bulkExprs), "_"), "[[", 1))
  
  cat("Dimension of psuedo-bulk expression")
  print(dim(bulkExprs))
  
  bulk_clustering_res <- lapply(aggExprs, function(x) {
    res <- unlist(lapply(strsplit(rownames(x), "_"), "[[", 2))
    res <- as.numeric(res)
    names(res) <- rownames(x)
    res
  })
  
  bulk_clustering_distProp <- lapply(bulk_clustering_res, function(x) {
    res <- rep(1, length(x))
    names(res) <- names(x)
    res
  })
  
  
  print("Identifying MNC using psuedo-bulk")
  
  
  if (sample_mode & match_within_batch) {
    print("match sample within batch")
    
    # sample_to_batch_idx <- match(sample_list[as.numeric(pseudo_batch)], meta_sample[, "sample"])
    # bulkExprs_batch <- meta_sample[sample_to_batch_idx, "batch"]
    bulk_clustering_res_batch <- meta_sample[match(pseudobulk_sample_list, meta_sample[, "sample"]), "batch"]
    
    replicate_vector_batch <- lapply(batch_list, function(x) {
      idx <- bulk_clustering_res_batch == x
      bulk_clustering_subset <- bulk_clustering_res[idx]
      mnc_res <- findMNC(bulkExprs[chosen.hvg, names(unlist(bulk_clustering_subset))], 
                         bulk_clustering_subset, dist = "cor",
                         BPPARAM = use_bpparam, 
                         plot_igraph = TRUE)
      
      replicate_vector <- mncReplicate(clustering_list = bulk_clustering_subset, 
                                       clustering_distProp = bulk_clustering_distProp[idx], 
                                       replicate_prop = 1, mnc_df = mnc_res)
      
      return(replicate_vector)
    })
    
    
    
    bulk_clustering_res_batch <- lapply(replicate_vector_batch, function(x) {
      res <- as.numeric(factor(x))
      names(res) <- names(x)
      res
    })
    
    
    mnc_res <- findMNC(bulkExprs[chosen.hvg, ], 
                       bulk_clustering_res_batch, dist = "cor",
                       BPPARAM = use_bpparam, 
                       plot_igraph = TRUE)
    
    
    bulk_clustering_distProp_batch <- lapply(bulk_clustering_res_batch, function(x) {
      res <- rep(1, length(x))
      names(res) <- names(x)
      res
    })
    
    
    
    
    replicate_vector <- mncReplicate(clustering_list = bulk_clustering_res_batch, 
                                     clustering_distProp = bulk_clustering_distProp_batch, 
                                     replicate_prop = 1, mnc_df = mnc_res)
    replicate_vector <- replicate_vector[colnames(bulkExprs)]
    
  }   
  
  if (sample_mode & !match_within_batch) {
    
    
    if (match_by_sample) {
      
      print("match by sample")
      
      mnc_res <- findMNC(bulkExprs[chosen.hvg, ], 
                         bulk_clustering_res, dist = "cor",
                         BPPARAM = use_bpparam, 
                         plot_igraph = TRUE)
      
      
      
      replicate_vector <- mncReplicate(clustering_list = bulk_clustering_res, 
                                       clustering_distProp = bulk_clustering_distProp, 
                                       replicate_prop = 1, mnc_df = mnc_res)
    } else {
      
      print("not match by sample")
      
      bulk_clustering_res_batch <- meta_sample[match(pseudobulk_sample_list, meta_sample[, "sample"]), "batch"]
      bulk_clustering_res_byBatch <- sapply(batch_list, function(x) {
        unlist(bulk_clustering_res[bulk_clustering_res_batch == x])
      })
      
      bulk_clustering_distProp_byBatch <- sapply(batch_list, function(x) {
        unlist(bulk_clustering_distProp[bulk_clustering_res_batch == x])
      })
      
      
      mnc_res <- findMNC(bulkExprs[chosen.hvg, ], 
                         bulk_clustering_res_byBatch, 
                         dist = "cor",
                         BPPARAM = use_bpparam, 
                         plot_igraph = TRUE)
      
      
      
      replicate_vector <- mncReplicate(clustering_list = bulk_clustering_res_byBatch, 
                                       clustering_distProp = bulk_clustering_distProp_byBatch, 
                                       replicate_prop = 1, mnc_df = mnc_res)
    }
    
    
    
  }
  
  if (cond_mode) {
    
    print("condition_mode")
    # sample_to_batch_idx <- match(sample_list[as.numeric(pseudo_batch)], meta_sample[, "sample"])
    # bulkExprs_batch <- meta_sample[sample_to_batch_idx, "batch"]
    bulk_clustering_res_condition <- meta_sample[match(pseudobulk_sample_list, meta_sample[, "sample"]), "condition"]
    
    replicate_vector_condition <- lapply(condition_list, function(x) {
      idx <- bulk_clustering_res_condition == x
      bulk_clustering_subset <- bulk_clustering_res[idx]
      mnc_res <- findMNC(bulkExprs[chosen.hvg, names(unlist(bulk_clustering_subset))], 
                         bulk_clustering_subset, dist = "cor",
                         BPPARAM = use_bpparam, 
                         plot_igraph = TRUE)
      
      replicate_vector <- mncReplicate(clustering_list = bulk_clustering_subset, 
                                       clustering_distProp = bulk_clustering_distProp[idx], 
                                       replicate_prop = 1, mnc_df = mnc_res)
      
    })
    
    replicate_vector_condition <- lapply(1:length(replicate_vector_condition), function(i) {
      res <- paste(condition_list[i], replicate_vector_condition[[i]], sep = "_")
      names(res) <- names(replicate_vector_condition[[i]])
      res
    })
    
    replicate_vector_condition <- unlist(replicate_vector_condition)
    replicate_vector <- replicate_vector_condition[colnames(bulkExprs)]
  }
  
  if (!cond_mode & !sample_mode) {
    
    mnc_res <- findMNC(bulkExprs[chosen.hvg, ], 
                       bulk_clustering_res, dist = "cor",
                       BPPARAM = use_bpparam, 
                       plot_igraph = TRUE)
    
    
    
    replicate_vector <- mncReplicate(clustering_list = bulk_clustering_res, 
                                     clustering_distProp = bulk_clustering_distProp, 
                                     replicate_prop = 1, mnc_df = mnc_res)
  }
  
  
  
  
  if (return_subset) {
    subset <- chosen.hvg
  } else {
    subset <- NULL
  }
  
  print("Running RUV")
  bulkExprs[is.infinite(bulkExprs)] <- 0
  bulkExprs <- as(bulkExprs, "dgCMatrix")
  res <- pseudoRUVIII(Y = t(exprsMat), 
                      Y_pseudo = t(bulkExprs), 
                      #batch = pseudo_batch,
                      M = ruv::replicate.matrix(replicate_vector), 
                      ctl = rownames(exprsMat) %in% ctl,
                      k = ruvK, 
                      BSPARAM = use_bsparam,
                      return.info = TRUE,
                      subset = subset,
                      normalised = normalised)
  
  
  
  
  gc(reset = TRUE)
  
  # return(t(res))
  return(res)
  
}






###################################################################################
##################          Pseudobulk                    #########################
###################################################################################


create_pseudoBulk_bagging <- function(exprs_old,
                                      reduced_coordinates,
                                      k = 50, 
                                      prop = 0.8) {
  
  # The following codes are from Cicero
  
  nn_map <- FNN::knn.index(reduced_coordinates, 
                           k=(k-1)) # no data.frame wrapper
  
  row.names(nn_map) <- row.names(reduced_coordinates)
  
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
  
  good_choices <- seq_len(nrow(nn_map))
  choice <- sample(seq_len(length(good_choices)), 
                   size = 1, replace = FALSE)
  chosen <- good_choices[choice]
  good_choices <- good_choices[good_choices != good_choices[choice]]
  
  it <- 0
  k2 <- k * 2 # Compute once
  
  # function for sapply
  get_shared <- function(other, this_choice) {
    k2 - length(union(cell_sample[other,], this_choice))
  }
  
  while (length(good_choices) > 0 & it < 1000) { # slow
    it <- it + 1
    choice <- sample(seq_len(length(good_choices)), 
                     size = 1, replace = FALSE)
    new_chosen <- c(chosen, good_choices[choice])
    good_choices <- good_choices[good_choices != good_choices[choice]]
    cell_sample <- nn_map[new_chosen,]
    
    others <- seq_len(nrow(cell_sample) - 1)
    this_choice <- cell_sample[nrow(cell_sample),]
    shared <- sapply(others, get_shared, this_choice = this_choice)
    
    if (max(shared) < prop * k) {
      chosen <- new_chosen
    }
  }
  
  cell_sample <- nn_map[chosen,]
  
  
  mask <- sapply(seq_len(nrow(cell_sample)), 
                 function(x) seq_len(ncol(exprs_old)) %in% 
                   cell_sample[x,,drop=FALSE])
  mask <- Matrix::Matrix(mask)
  mask <- mask/nrow(cell_sample)
  new_exprs <- exprs_old %*% mask
  
  print(dim(new_exprs))
  
  return(new_exprs)
  
  
}





create_pseudoBulk <- function(exprsMat, cell_info, k_fold = 20, ncore = 8) {
  
  k_fold <- min(ncol(exprsMat), k_fold)
  cv <- cvTools::cvFolds(ncol(exprsMat), K = k_fold)
  
  exprsMat_pseudo <- parallel::mclapply(seq_len(k_fold), function(i) {
    subset_idx <- cv$subsets[cv$which == i]
    cellType_tab <- table(droplevels(factor(cell_info[subset_idx])))
    cellTypes_n_mat <- matrix(rep(cellType_tab, 
                                  nrow(exprsMat)), 
                              nrow = length(cellType_tab), byrow = FALSE)
    
    
    system.time(res<- Matrix.utils::aggregate.Matrix(t(exprsMat[, subset_idx]), 
                                                     cell_info[subset_idx],
                                                     fun = "sum"))
    res <- res/cellTypes_n_mat
    
    rownames(res) <- paste(rownames(res),
                           i, sep = "_")
    res
  }, mc.cores = ncore)
  
  exprsMat_pseudo <- do.call(rbind, exprsMat_pseudo)
  
  return(exprsMat_pseudo)
}





create_pseudoBulk_sum <- function(exprsMat, cell_info, k_fold = 20, ncore = 8) {
  
  k_fold <- min(ncol(exprsMat), k_fold)
  cv <- cvTools::cvFolds(ncol(exprsMat), K = k_fold)
  
  exprsMat_pseudo <- parallel::mclapply(seq_len(k_fold), function(i) {
    subset_idx <- cv$subsets[cv$which == i]
    cellType_tab <- table(droplevels(factor(cell_info[subset_idx])))
    cellTypes_n_mat <- matrix(rep(cellType_tab, 
                                  nrow(exprsMat)), 
                              nrow = length(cellType_tab), byrow = FALSE)
    
    
    system.time(res<- Matrix.utils::aggregate.Matrix(t(exprsMat[, subset_idx]), 
                                                     cell_info[subset_idx],
                                                     fun = "sum"))
    res <- scater::normalizeCounts(res)
    
    rownames(res) <- paste(rownames(res),
                           i, sep = "_")
    res
  }, mc.cores = ncore)
  
  exprsMat_pseudo <- do.call(rbind, exprsMat_pseudo)
  
  return(exprsMat_pseudo)
}




create_pseudoBulk_sample <- function(exprsMat, cell_info, prop_pseudoBulk = 0.1, 
                                     ncore = 8,
                                     upsample = FALSE,
                                     min_cells) {
  
  if (upsample) {
    prop_pseudoBulk <- max(prop_pseudoBulk, min_cells/ncol(exprsMat))
    min_cells_clust <- min(table(cell_info))
  }
  
  
  exprsMat_pseudo <- lapply(split(seq_len(ncol(exprsMat)), cell_info), function(x) {
    if (length(x) < 5) {
      num_cells <- min(length(x), 5)
    } else {
      if (upsample) {
        num_cells <- max(min_cells_clust, round(prop_pseudoBulk * length(x)))
      } else {
        num_cells <- round(prop_pseudoBulk * length(x))
      }
      num_cells <- max(num_cells, 5)
    }
    
    tmp <- exprsMat[, sample(x, num_cells), drop = FALSE]
    colnames(tmp) <- seq_len(ncol(tmp))
    tmp
  })
  
  exprsMat_pseudo <- lapply(1:length(exprsMat_pseudo), function(i) {
    colnames(exprsMat_pseudo[[i]]) <- paste(i, colnames(exprsMat_pseudo[[i]]), sep = "_")
    exprsMat_pseudo[[i]]
  })
  
  
  exprsMat_pseudo <- do.call(cbind, exprsMat_pseudo)
  
  return(t(exprsMat_pseudo))
}



create_pseudoBulk_kmeans <- function(exprsMat, cell_info, hvg = NULL, k_fold = 20, ncore = 8) {
  
  if (is.null(hvg)) {
    hvg <- rownames(exprsMat)
  }
  
  k_fold <- min(ncol(exprsMat), k_fold)
  
  
  exprsMat_pseudo <- lapply(split(seq_len(ncol(exprsMat)), cell_info), function(x) {
    pca <- irlba::prcomp_irlba(t(exprsMat[hvg, x]), n = 10)
    km <- kmeans(pca$x, centers = k_fold, nstart = 10, iter.max = 50)
    
    cellType_tab <- table(km$cluster)
    
    cellTypes_n_mat <- matrix(rep(cellType_tab, 
                                  ncol(t(exprsMat[, x]))), 
                              nrow = length(cellType_tab), byrow = FALSE)
    res<- Matrix.utils::aggregate.Matrix(t(exprsMat[, x]), 
                                         km$cluster,
                                         fun = "sum")
    res <- t(res/cellTypes_n_mat)
    
  })
  
  exprsMat_pseudo <- lapply(1:length(exprsMat_pseudo), function(i) {
    colnames(exprsMat_pseudo[[i]]) <- paste(i, colnames(exprsMat_pseudo[[i]]), sep = "_")
    exprsMat_pseudo[[i]]
  })
  
  
  exprsMat_pseudo <- do.call(cbind, exprsMat_pseudo)
  
  return(t(exprsMat_pseudo))
}





create_pseudoBulk_graph <- function(exprsMat, cell_info, hvg = NULL, k_fold = 20, ncore = 2) {
  
  if (is.null(hvg)) {
    hvg <- rownames(exprsMat)
  }
  
  cell_info <- as.character(cell_info)
  exprsMat_pseudo <- parallel::mclapply(split(seq_len(ncol(exprsMat)), cell_info), function(x) 
  {
    if (length(x) > 5) {
      g <- suppressWarnings(buildKNNGraph(exprsMat[, x], 
                                          k = min(k_fold, round(length(x)/2)), 
                                          BNPARAM = BiocNeighbors::AnnoyParam(),
                                          subset.row = hvg))
      cluster_res <- igraph::cluster_louvain(g)
      cellType_tab <- table(cluster_res$membership)
      cellTypes_n_mat <- matrix(rep(cellType_tab, 
                                    ncol(t(exprsMat[, x]))), 
                                nrow = length(cellType_tab), byrow = FALSE)
      res <- Matrix.utils::aggregate.Matrix(t(exprsMat[, x]), 
                                            cluster_res$membership,
                                            fun = "sum")
      res <- t(res/cellTypes_n_mat)
    } else {
      res <- exprsMat[, x, drop = FALSE]
      colnames(res) <- seq_len(ncol(res))
      res
    }
    
    
  }, mc.cores = ncore)
  #print(lapply(exprsMat_pseudo, head))
  exprsMat_pseudo <- lapply(1:length(exprsMat_pseudo), function(i) {
    colnames(exprsMat_pseudo[[i]]) <- paste(i, colnames(exprsMat_pseudo[[i]]), sep = "_")
    exprsMat_pseudo[[i]]
  })
  
  
  exprsMat_pseudo <- do.call(cbind, exprsMat_pseudo)
  
  return(t(exprsMat_pseudo))
}




create_pseudoBulk_pool <- function(exprsMat_counts, cell_info, k_fold = 10, ncore = 2) {
  
  
  cell_info <- as.character(cell_info)
  exprsMat_pseudo <- parallel::mclapply(split(seq_len(ncol(exprsMat_counts)), cell_info), function(x) 
  {
    if (length(x) > k_fold) {
      # pool
      total_counts <- colSums(exprsMat_counts[, x])
      bins <- cut(rank(total_counts), k_fold)
      res <- Matrix.utils::aggregate.Matrix(t(exprsMat_counts[, x, drop = FALSE]), 
                                            bins,
                                            fun = "sum")
      bins_tab <- table(bins)
      cellTypes_n_mat <- matrix(rep(bins_tab[rownames(res)], 
                                    ncol(t(exprsMat_counts[, x]))), 
                                nrow = nrow(res), byrow = FALSE)
      res <- t(res/cellTypes_n_mat)
      res <- scuttle::normalizeCounts((res), log = TRUE)
      
    } else {
      res <- exprsMat_counts[, x, drop = FALSE]
      res <- scuttle::normalizeCounts((res), log = TRUE)
      colnames(res) <- seq_len(ncol(res))
      res
    }
    
    
  }, mc.cores = ncore)
  #print(lapply(exprsMat_pseudo, head))
  exprsMat_pseudo <- lapply(1:length(exprsMat_pseudo), function(i) {
    colnames(exprsMat_pseudo[[i]]) <- paste(i, colnames(exprsMat_pseudo[[i]]), sep = "_")
    exprsMat_pseudo[[i]]
  })
  
  
  exprsMat_pseudo <- do.call(cbind, exprsMat_pseudo)
  
  return(t(exprsMat_pseudo))
}





create_pseudoBulk_pool_divide <- function(exprsMat_counts, cell_info, k_fold = 10, ncore = 2) {
  
  
  cell_info <- as.character(cell_info)
  exprsMat_pseudo <- parallel::mclapply(split(seq_len(ncol(exprsMat_counts)), cell_info), function(x) 
  {
    if (length(x) > k_fold) {
      # pool
      total_counts <- colSums(exprsMat_counts[, x])
      bins <- cut(rank(total_counts), k_fold)
      res <- Matrix.utils::aggregate.Matrix(t(exprsMat_counts[, x, drop = FALSE]), 
                                            bins,
                                            fun = "sum")
      bins_tab <- table(bins)
      cellTypes_n_mat <- matrix(rep(bins_tab[rownames(res)], 
                                    ncol(t(exprsMat_counts[, x]))), 
                                nrow = nrow(res), byrow = FALSE)
      res <- lapply(seq_len(nrow(res)), function(k) rbinom(ncol(res), round(res[k, ]), 1/cellTypes_n_mat[k, ]))
      res <- do.call(cbind, res)
      res <- scater::normalizeCounts((res))
      colnames(res) <- seq_len(ncol(res))
      res
    } else {
      res <- exprsMat_counts[, x, drop = FALSE]
      res <- scuttle::normalizeCounts((res), log = TRUE)
      colnames(res) <- seq_len(ncol(res))
      res
    }
    
    
  }, mc.cores = ncore)
  
  exprsMat_pseudo <- lapply(1:length(exprsMat_pseudo), function(i) {
    colnames(exprsMat_pseudo[[i]]) <- paste(i, colnames(exprsMat_pseudo[[i]]), sep = "_")
    exprsMat_pseudo[[i]]
  })
  
  
  exprsMat_pseudo <- do.call(cbind, exprsMat_pseudo)
  
  return(t(exprsMat_pseudo))
}



create_pseudoBulk_divide <- function(exprsMat_counts, cell_info, k_fold = 20, ncore = 8) {
  
  
  
  
  exprsMat_pseudo <- parallel::mclapply(split(seq_len(ncol(exprsMat_counts)), cell_info), function(x) 
  {
    if (length(x) > k_fold) {
      # pool
      
      cv <- cvTools::cvFolds(length(x), K = k_fold)
      bins <- cv$which[cv$subsets]
      res <- Matrix.utils::aggregate.Matrix(t(exprsMat_counts[, x, drop = FALSE]), 
                                            bins,
                                            fun = "sum")
      bins_tab <- table(bins)
      cellTypes_n_mat <- matrix(rep(bins_tab, 
                                    ncol(t(exprsMat_counts[, x]))), 
                                nrow = length(bins_tab), byrow = FALSE)
      res <- lapply(seq_len(nrow(res)), function(k) rbinom(ncol(res),  round(res[k, ]), 1/cellTypes_n_mat[k, ]))
      res <- do.call(cbind, res)
      res <- scater::normalizeCounts((res))
      colnames(res) <- as.numeric(factor(names(bins_tab)))
      res
    } else {
      res <- exprsMat_counts[, x, drop = FALSE]
      colnames(res) <- seq_len(ncol(res))
      res
    }
    
    
  }, mc.cores = ncore)
  
  
  exprsMat_pseudo <- lapply(1:length(exprsMat_pseudo), function(i) {
    colnames(exprsMat_pseudo[[i]]) <- paste(i, colnames(exprsMat_pseudo[[i]]), sep = "_")
    exprsMat_pseudo[[i]]
  })
  
  
  exprsMat_pseudo <- do.call(cbind, exprsMat_pseudo)
  
  
  return(t(exprsMat_pseudo))
}




create_pseudoBulk_graph_divide <- function(exprsMat, exprsMat_counts, cell_info, hvg = NULL, k_fold = 20, ncore = 2) {
  
  if (is.null(hvg)) {
    hvg <- rownames(exprsMat)
  }
  
  cell_info <- as.character(cell_info)
  exprsMat_pseudo <- parallel::mclapply(split(seq_len(ncol(exprsMat)), cell_info), function(x) 
  {
    if (length(x) > 5) {
      g <- suppressWarnings(buildKNNGraph(exprsMat[, x], 
                                          k = min(k_fold, round(length(x)/2)), 
                                          BNPARAM = BiocNeighbors::AnnoyParam(),
                                          subset.row = hvg))
      cluster_res <- igraph::cluster_louvain(g)$membership
      
      total_counts <- colSums(exprsMat_counts[, x])
      res <- Matrix.utils::aggregate.Matrix(t(exprsMat_counts[, x, drop = FALSE]), 
                                            cluster_res,
                                            fun = "sum")
      cluster_tab <- table(cluster_res)
      cellTypes_n_mat <- matrix(rep(cluster_tab, 
                                    ncol(t(exprsMat_counts[, x]))), 
                                nrow = nrow(res), byrow = FALSE)
      res <- lapply(seq_len(nrow(res)), function(k) rbinom(ncol(res), res[k, ], 1/cellTypes_n_mat[k, ]))
      res <- do.call(cbind, res)
      res <- scater::normalizeCounts((res))
      colnames(res) <- as.numeric(factor(names(cluster_tab)))
      res
    } else {
      res <- exprsMat[, x, drop = FALSE]
      colnames(res) <- seq_len(ncol(res))
      res
    }
    
    
  }, mc.cores = ncore)
  #print(lapply(exprsMat_pseudo, head))
  exprsMat_pseudo <- lapply(1:length(exprsMat_pseudo), function(i) {
    colnames(exprsMat_pseudo[[i]]) <- paste(i, colnames(exprsMat_pseudo[[i]]), sep = "_")
    exprsMat_pseudo[[i]]
  })
  
  
  exprsMat_pseudo <- do.call(cbind, exprsMat_pseudo)
  
  return(t(exprsMat_pseudo))
}




create_pseudoBulk_KNN <- function(x, cell_info, prop_pseudoBulk = 0.1, 
                                  chosen.hvg, pseudobulk_k = 10,
                                  upsample = FALSE,
                                  min_cells = NULL) {
  
  if (is.null(pseudobulk_k)) {
    pseudobulk_k <- 20
  }
  
  if (ncol(x) < 50) {
    
    exprsMat_pseudo <- t(x)
    rownames(exprsMat_pseudo) <- paste(cell_info, seq_len(nrow(exprsMat_pseudo)), sep = "_")
    return(exprsMat_pseudo)
    
  } else {
    if (upsample) {
      num_cells <-  max(min_cells, round(prop_pseudoBulk * ncol(x)))
    } else {
      num_cells <-  round(prop_pseudoBulk * ncol(x))
      
    }
    
    
    out <- BiocNeighbors::findKNN(t(x[chosen.hvg, ]), k = pseudobulk_k,
                                  BNPARAM = BiocNeighbors::AnnoyParam())
    idx <- sample(ncol(x), num_cells)
    out <- reshape2::melt((out$index[idx, ]))
    
    mat <- x[, out$value]
    
    cellType_tab <- table(droplevels(factor(out$Var1)))
    cellTypes_n_mat <- matrix(rep(cellType_tab, 
                                  nrow(mat)), 
                              nrow = length(cellType_tab), byrow = FALSE)
    
    exprsMat_pseudo <- Matrix.utils::aggregate.Matrix(t(mat), 
                                                      out$Var1,
                                                      fun = "sum")
    exprsMat_pseudo <- exprsMat_pseudo/cellTypes_n_mat
    
    
    # exprsMat_pseudo <- t(exprsMat_pseudo)
    # 
    rownames(exprsMat_pseudo) <- paste(cell_info[idx], rownames(exprsMat_pseudo), sep = "_")
    return(exprsMat_pseudo)
  }
  
  
}

