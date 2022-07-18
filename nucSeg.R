
nucSeg <- function(image,
                   nucleus_index = 1,
                   size_selection = 10,
                   smooth = 1,
                   tolerance = 0.01,
                   watershed = "combine",
                   ext = 1) {
  
  
  
  
  # Prepare matrix use to segment nuclei
  nuc <- .prepNucSignal(image, nucleus_index, smooth)
  
  # Segment Nuclei
  nth <-
    EBImage::otsu(nuc, range = range(nuc))  # thresholding on the sqrt intensities works better.
  nMask <-
    nuc > nth  # the threshold is squared to adjust of the sqrt previously.
  
  
  # Size Selection
  nMaskLabel <- EBImage::bwlabel(nMask*1)
  tabNuc <- table(nMaskLabel)
  nMask[nMaskLabel %in% names(which(tabNuc <= size_selection))] <- 0  
  nMaskLabel[nMaskLabel %in% names(which(tabNuc <= size_selection))] <- 0 
  
  
  if(watershed == "distance"){
    if(is.null(tolerance)) tolerance <- 1
    dist <- EBImage::distmap(nMask)
    wMask <-
      EBImage::watershed(dist, tolerance = tolerance, ext = ext)
      return(wMask)
    
  }
  
  if(watershed == "combine"){
    
    # Scale cell intensities
    avg <- tapply(nuc, nMaskLabel, mean)
    AVG <- nMask
    AVG[] <- avg[as.character(nMaskLabel)]
    nuc <- (nuc/AVG)*nMask
    nuc <- nuc/median(avg[as.character(nMaskLabel)])*nMask
    #nuc <- nuc/sd(nuc[nuc!=0])/2
    
    dist <- EBImage::distmap(nMask)
    
    if(is.null(tolerance)){
      tolerance <- .estimateTolerance(dist*nuc, nMask)
    }
    
    wMask <-
      EBImage::watershed(dist*nuc, tolerance = tolerance, ext = ext)
    return(wMask)
  }
  
  
  # Add distance to nuc signal
  cellRadius <- 2*floor(sqrt(size_selection/pi)/2)+1
  nuc <- filter2(nuc, makeBrush(cellRadius, shape='disc'))
  nuc <- nuc * nMask
  
  if(is.null(tolerance)){
    tolerance <- .estimateTolerance(nuc, nMask)
  }
  
  wMask <-
    EBImage::watershed(nuc, tolerance = tolerance, ext = ext)
  
  # Size Selection
  tabNuc <- table(wMask)
  wMask[wMask %in% names(which(tabNuc <= size_selection))] <- 0  

  wMask
  
}


.prepNucSignal <- function(image, nucleus_index, smooth){
  
  
  if("PCA" %in% nucleus_index){
  image <- apply(image, 3, function(x){
    x <- (x)
    EBImage::gblur(x, smooth)
  }, simplify = FALSE)
  
  image <- abind(image, along = 3)
  
  image.long <- apply(image,3, as.numeric)
  pca <- prcomp((image.long))
  
  usePC <- 1
  if(any(nucleus_index%in%colnames(image.long))){
    ind <- intersect(nucleus_index, colnames(image.long))
    usePC <- which.max(apply(pca$x, 2, cor, image.long[,nucleus_index[nucleus_index != "PCA"][1]]))
  }
  
  imagePC <- as.matrix(image[,,1])
  imagePC[] <- pca$x[,usePC] - min(pca$x[,usePC])
  return(imagePC)
  }
  
  if(is(nucleus_index, "character"))
  ind <- intersect(nucleus_index, dimnames(image)[[3]])
  
  nuc <- image[, , ind]
  if(length(ind)>1) nuc <- apply(nuc, c(1,2), mean)
  nuc <- EBImage::gblur(nuc, smooth)
  
  nuc
  
}


.estimateTolerance <- function(input, nMask){
  y <- EBImage::distmap(nMask)
   fit <- lm(as.numeric(input[y>0]) ~ as.numeric(y[y>0])-1)
   tolerance <- coef(fit)[1]
   tolerance
  # tolerance <- sd(as.numeric(input[y>0]))/sd(as.numeric(y[y>0]))
  # tolerance
}
 
