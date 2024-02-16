#' Apply Mosaic-PICASSO Correction to an Image
#'
#' This function orchestrates the entire Mosaic-PICASSO algorithm to infer and remove autofluorescence
#' and spillover from a multiplex image. Steps include loading the image, splitting the image into tiles,
#' selecting only tiles with the lowest structural similarity index measure, calculating the matrix of scaling
#' factors for correction, applying the correction to each channel, optionally enhancing the contrast of
#' the corrected image, and saving the output. The input image must be an 8-bit grayscale multichannel .tif.
#' The output is an image of the same format, with "_NEW" appended to its name, along with a .csv copy of
#' the correction matrix.
#'
#' @description Apply Mosaic-PICASSO Correction to an Image
#' @param file        Character: path to 8-bit grayscale multichannel .tif image
#' @param method      Character: whether to minimize mutual info "MI" or correlation "Cor"
#' @param q_thr       Numeric: quantile threshold to select only the tiles with lowest SSIM
#' @param tile_size   Integer: dimensions of each square tile in pixels
#' @param stride      Integer: distance between top left corner of each tile in pixels
#' @param pix_ss      Integer: number of pixels to randomly subsample for calculating MI or correlation
#' @param enhance     Logical: whether to enhance contrast in the corrected image
#' @param groups      List: vectors of channel numbers to include in self-contained groups for correction
#' @param cycles      Integer: number of iterations of correction to perform
#' @param cores       Integer: number of cores to use for parallel calculations
#' @return      Nothing: a corrected .tif and a .csv of the correction matrix are saved to the working directory
#' @export
mosaic_picasso <- function(file, method="MI", q_thr=0.5, tile_size=50, stride=50, pix_ss="all", enhance=T, 
                           groups=NULL, cycles=1, cores=NULL) {
  print("Loading image.")
  IMG <- load_image(file)                                                       # Load img to 3d array format
  N_CH <- dim(IMG)[3]                                                           # Number of channels
  if (enhance) {                                                                # If image will be contrast-enhanced
    frc_sat <- rep(NA, N_CH)                                                    # Get fraction of saturated pix
    for (i in 1:N_CH) {                                                         # in each channel
      frc_sat[i] <- sum(IMG[,,i] == 255) / prod(dim(IMG)[1:2])
    }
  }
  if (!is.null(groups)) {
    print("Splitting image into channel groups for correction.")
    if (any(table(unlist(groups)) > 1)) {
      stop("Error: Each channel can only be used in one group.")
    }
    if (length(unlist(groups)) != N_CH) {
      print("Warning: Channels that are not included in a group will not be used for correction.")
    }
    img <- vector("list", length(groups))
    mat <- vector("list", length(groups))
    names(mat) <- paste("G",1:length(groups))
    for (g in 1:length(groups)) {
      img[[g]] <- IMG[,,groups[[g]],drop=F]
    }
  } else {
    img <- vector("list", 1)
    img[[1]] <- IMG
    mat <- vector("list", 1)
  }
  for (g in 1:length(img)) {
    print(paste("Beginning channel group ", g, ".", sep=""))
    n_ch <- dim(img[[g]])[3]                                                    # Number of channels in this group
    p_mat_0 <- diag(n_ch)                                                       # Initial P matrix                                                        
    p_mat_1 <- diag(n_ch)                                                       # Updated P matrix
    print("Calculating correction matrix.")
    for (i in 1:cycles) {                                                       # For each cycle of correction
      print(paste("   Cycle ", i,".",sep=""))
      if (i > 1) {                                                              # Correct image via initial p matrix
        tmp_img <- correct_image(img[[g]], p_mat_0)
      } else {                                                                  # Skip this step if this is 1st cycle
        tmp_img <- img[[g]]
      }
      tmp_mat <- calc_Pmat(tmp_img, method, q_thr, tile_size, stride, pix_ss, cores)# New P matrix from corrected image 
      p_mat_1 <- update_Pmat(p_mat_0, tmp_mat)                                  # Updated P matrix = initial * new
      p_mat_0 <- p_mat_1                                                        # Updated P matrix becomes initial
    }
    print("Applying correction matrix.")
    img[[g]] <- correct_image(img[[g]], p_mat_1)                                # Correct image via updated P matrix
    mat[[g]] <- p_mat_1
  }
  for (i in 1:N_CH) {
    grp <- which(unlist(lapply(groups, function(x,y) {y %in% x}, y=i)))
    idx <- which(groups[[grp]] == i)
    IMG[,,i] <- img[[grp]][,,idx]
  }
  if (enhance) {                                                                # If contrast-enhancing final image
    print("Enhancing contrast.")
    for (i in 1:N_CH) {
      if (frc_sat[i] == 0) {
        mult <- 255
      } else {
        mult <- unname(quantile(c(IMG[,,i]), 1-frc_sat[i]))                     # Multiplier to brighten image
      }
      if (mult > 0) {                                                           # If the image is not totally black
        IMG[,,i] <- round(IMG[,,i] * (255/mult))                                # Brighten image and round pixels
        IMG[,,i][IMG[,,i] > 255] <- 255                                         # Any pixels > 255 are set to 255
      }
    }
  }
  print("Saving corrected image.")
  out_file <- paste(sub(".tif", "", file), "_NEW.tif", sep="")                  # Add _NEW to file name
  ijtiff::write_tif(IMG, out_file, overwrite=T)                                 # Save corrected image
  print("Saving correction matrix.")
  out_file <- paste(sub(".tif", "", file), "_MAT.xlsx", sep="")                  # Add _MAT to file name
  openxlsx::write.xlsx(mat, out_file)                                           # Save correction matrix
  print("Finished.")
}

#' Load an Image
#'
#' This function loads an 8-bit grayscale multichannel .tif image and transforms it into a 3d array.
#'
#' @description Load an image
#' @param file        Character: path to 8-bit grayscale multichannel .tif image
#' @return      Array 3d: a 3d array representation of the original image
load_image <- function(file) {
  img <- tiff::readTIFF(file, all = T, as.is = F)                               # Initial read-in with tiff package
  img_dim <- c(dim(img[[1]])[1], dim(img[[1]])[2], length(img))                 # Image is list of channels
  for (i in 1:length(img)) {                                                    # For each channel
    img[[i]] <- img[[i]]*255                                                    # get pixel values in 0-255
    mode(img[[i]]) <- "integer"                                                 # and store in integer format
  }
  img <- simplify2array(img)                                                    # Collapse channel list into 3d array
  return(img)
} 

#' Assemble the Correction Matrix
#'
#' This function orchestrates the inference of the correction scaling factor for every pair of channels.
#' Entry i,j in the correction matrix gives the scaling factor to apply to image i before subtracting it
#' from image j, in order to correct image j based on overlap with image i.
#'
#' @description Assemble the correction matrix
#' @param img         Array 3d: all channels of an image
#' @param method      Character: whether to minimize mutual info "MI" or correlation "Cor"
#' @param q_thr       Numeric: quantile threshold to select only the tiles with lowest SSIM
#' @param tile_size   Integer: dimensions of each square tile in pixels
#' @param stride      Integer: distance between top left corner of each tile in pixels
#' @param pix_ss      Integer: number of pixels to randomly subsample for calculating MI or correlation
#' @param cores       Integer: number of cores to use for parallel calculations
#' @return      Matrix: the fully calculated correction matrix
calc_Pmat <- function(img, method, q_thr, tile_size, stride, pix_ss, cores) {
  n_ch <- dim(img)[3]
  p_mat <- diag(n_ch)
  idx <- expand.grid(1:n_ch, 1:n_ch)
  idx <- idx[idx[,1] != idx[,2], ]
  img_chp <- chop_image(img, tile_size, stride)
  num_cores <- min(parallel::detectCores() - 1, cores)                    
  make_cluster <- parallel::makeCluster(num_cores, type = "PSOCK")    
  doParallel::registerDoParallel(cl = make_cluster)            
  `%dopar%` <- foreach::`%dopar%`                              
  a <- foreach::foreach (i = c(1:nrow(idx)), .export=c('calc_scalar','calc_ssim','calc_MI','calc_corr')) %dopar% {
    calc_scalar(img_chp, idx[i,1], idx[i,2], method, q_thr, pix_ss)
  }
  parallel::stopCluster(cl = make_cluster) 
  a <- unlist(a)
  k <- 1
  for (i in 1:nrow(idx)) {
    if (idx[i,1] != idx[i,2]) {
      p_mat[idx[i,1], idx[i,2]] <- -a[k]
      k <- k + 1
    }
  }
  return(p_mat)
}

#' Convert an Image into Tiles
#'
#' This function cuts an image into tiles of specified size and spacing. Each channel becomes a 3d array,
#' where x and y correspond to x and y in the original image, and z is the tile ID. The 3d arrays representing
#' each channel are returned in a list. 
#'
#' @description Convert an image into tiles
#' @param img         Array 3d: all channels of an image
#' @param tile_size   Integer: dimensions of each square tile in pixels
#' @param stride      Integer: distance between top left corner of each tile in pixels
#' @return      List: for each channel of the image, a 3d array of the tiles
chop_image <- function(img, tile_size, stride) {
  coors <- expand.grid(row = seq(0,floor(dim(img)[1]/stride)-1)*stride + 1,     # Coors of top-left
                       col = seq(0,floor(dim(img)[2]/stride)-1)*stride + 1)     # corner of each tile
  blocks <- kronecker(matrix(1:nrow(coors), nrow=length(unique(coors$row)),     # Nearly full img with
                             ncol=length(unique(coors$col)), byrow=F),          # all pix labeled by
                      matrix(1, stride, stride))                                # ID of their tile
  d1 <- dim(blocks)[1]                                                          # Width of block array
  d2 <- dim(blocks)[2]                                                          # Height of block array
  tiles <- vector("list", 2)                                                    # 2 sets of tiles
  for (i in 1:dim(img)[3]) {                                                    # For each channel
    tiles[[i]] <- lapply(split(img[1:d1,1:d2,i], blocks), matrix, nrow = stride)# List of 2d tiles
    tiles[[i]] <- simplify2array(tiles[[i]])                                    # 3d array of tiles
    if (tile_size < stride) {                                                   # If tiles should be
      tiles[[i]] <- tiles[[i]][1:tile_size, 1:tile_size, ]                      # smaller, cut them
    }
  }
  return(tiles)                       
}

#' Calculate the Scaling Factor for a Pair of Channels
#'
#' This function calculates the correction scaling factor for a single specified pair of channels.
#' The scaling factor will be applied to channel i before subtracting it from channel j, in order to 
#' correct channel j based on overlap with channel i. Candidate scaling factors are 0, 0.01, 0.02, ..., 1.00.
#' The chosen scaling factor A is the one that minimizes MI or |correlation| between channel i and 
#' channel j - A*Channel i. Only tiles below the specified quantile of SSIM between the two channels 
#' are considered when calculating MI or |correlation|. Among these tiles, even fewer pixels may be subsampled.
#'
#' @description Calculate the scaling factor for a pair of channels
#' @param img_chp     List: 3d array of tiles from chopped image, across all channels
#' @param i           Integer: index of 1st channel, which the 2nd channel is corrected with respect to
#' @param j           Integer: index of 2nd channel, which is used to correct the 1st image
#' @param method      Character: whether to minimize mutual info "MI" or correlation "Cor"
#' @param q_thr       Numeric: quantile threshold to select only the tiles with lowest SSIM
#' @param pix_ss      Integer: number of pixels to randomly subsample for calculating MI or correlation
#' @return      Numeric: the scaling factor to apply to channel i before subtracting from channel j
calc_scalar <- function(img_chp, i, j, method, q_thr, pix_ss) {
  img_chp <- img_chp[c(i,j)]
  ssim <- calc_ssim(img_chp)
  thr <- quantile(ssim, q_thr)
  img_chp_sub <- vector("list", 2)
  for (i in 1:2) {
    img_chp_sub[[i]] <- img_chp[[i]][,,ssim <= thr]
  }
  if (method == "MI") {
    obj_func <- function(a, img1, img2) {
      return(calc_MI(img1, img2 - a*img1, pix_ss))
    }
  } else if (method == "Cor") {
    obj_func <- function(a, img1, img2) {
      return(calc_corr(img1, img2 - a*img1, pix_ss))
    }
  }
  a_try <- as.list(seq(0,1,by=0.01))
  obj_out <- lapply(a_try, obj_func, img1=img_chp_sub[[1]], img2=img_chp_sub[[2]])
  out <- unlist(a_try)[which.min(unlist(obj_out))]
  if (length(out)==0) {
    return(0) # If cor/MI is NA for every alpha value, choose 0 as the alpha value
  } else {
    return(out)
  }
}

#' Calculate SSIM between Tiles of Two Channels
#'
#' This function calculates the Structural Similarity Index Measure, SSIM, between corresponding tiles
#' across two channels of an image. 
#'
#' @description Calculate SSIM between tiles of two channels
#' @param tiles   List: 3d arrays of tiles for two channels
#' @return      Vector: SSIM between corresponding pairs of tiles
calc_ssim <- function(tiles) {
  ss <- rep(NA, dim(tiles[[1]])[3])
  for (i in 1:length(ss)) {
    ss[i] <- SpatialPack::SSIM(tiles[[1]][,,i], tiles[[2]][,,i])$SSIM
  }
  return(ss)
}

#' Calculate Correlation of Pixel Values between Two Channels
#'
#' This function takes two channels, already restricted to just the correct subset of tiles, and 
#' calculates the pixel-wise correlation between the two channels. If pixels will be randomly subsampled,
#' this is done first. The absolute value of the correlation is returned.
#'
#' @description Calculate correlation of pixel values between two channels.
#' @param img1     Array 3d: tiles of the 1st channel for which correlation is sought with the 2nd
#' @param img2     Array 3d: tiles of the 2nd channel for which correlation is sought with the 1st
#' @param pix_ss   Integer: number of pixels to randomly subsample for calculating correlation
#' @return      Numeric: absolute value of correlation between the two channels
calc_corr <- function(img1, img2, pix_ss) {
  img1 <- as.vector(img1)
  img2 <- as.vector(img2)
  if (!is.character(pix_ss)) {
    pix_ss <- min(length(img1), pix_ss)
    keep_pix <- sample(length(img1), pix_ss)
    img1 <- img1[keep_pix]
    img2 <- img2[keep_pix]
  }
  if (length(unique(img1)) == 0 || length(unique(img2))==0) {
    return(NA) # If there is no variation in an img, cannot compute correlation.
  } else {
    cor_coef <- abs(cor(img1, img2))
    return(cor_coef)
  }
}

#' Calculate Mutual Information of Pixel Values between Two Channels
#'
#' This function takes two channels, already restricted to just the correct subset of tiles, and 
#' calculates the pixel-wise mutual information between the two channels. If pixels will be randomly subsampled,
#' this is done first. The normalized mutual information is returned.
#'
#' @description Calculate mutual information of pixel values between two channels.
#' @param img1     Array 3d: tiles of the 1st channel for which MI is sought with the 2nd
#' @param img2     Array 3d: tiles of the 2nd channel for which MI is sought with the 1st
#' @param pix_ss   Integer: number of pixels to randomly subsample for calculating MI
#' @return      Numeric: normalized mutual information between the two channels
calc_MI <- function(img1, img2, pix_ss) {
  img1 <- as.vector(img1)
  img2 <- as.vector(img2)
  if (!is.character(pix_ss)) {
    pix_ss <- min(length(img1), pix_ss)
    keep_pix <- sample(length(img1), pix_ss)
    img1 <- img1[keep_pix]
    img2 <- img2[keep_pix]
  }
  if (length(unique(img1)) == 0 || length(unique(img2))==0) {
    return(NA) # If there is no variation in an img, cannot compute mutual information.
  } else {
    img1 <- round(img1)
    img2 <- round(img2)
    img2 <- img2 + min(img2)
    ent1 <- unname(table(img1))
    ent2 <- unname(table(img2))
    ent1 <- ent1 / sum(ent1)
    ent2 <- ent2 / sum(ent2)
    ent1 <- -sum(ent1 * log2(ent1))
    ent2 <- -sum(ent2 * log2(ent2))
    imgb <- img1*(max(img2)+1) + img2
    entb <- unname(table(imgb)) 
    entb <- entb / sum(entb)
    entb <- -sum(entb * log2(entb))
    mi <- ent1 + ent2 - entb
    return((2*mi)/(ent1 + ent2))
  }
}

#' Update a Correction Matrix
#'
#' This function multiplies an initial correction matrix by a second one. This is only meaningful if
#' more than one cycle of correction is being performed.
#'
#' @description Update a correction matrix by multiplying an original correction matrix with a new one.
#' @param p_mat_0     Matrix: initial correction matrix
#' @param p_mat_1     Matrix: new correction matrix
#' @return      Matrix: updated correction matrix
update_Pmat <- function(p_mat_0, p_mat_1) {
  return(p_mat_1 %*% p_mat_0)
}

#' Apply Correction Matrix to an Image
#'
#' This function corrects a channel using a full correction matrix. Each channel is set to its original self, 
#' minus a scaled version of each of the other channels. With respect to the correction matrix, 
#' the column is the focal channel to be corrected, and the rows are the scaled factors for each of the other
#' channels to be subtracted.
#'
#' @description Apply the correction matrix of scaling factors to every pair of channels to correct an image.
#' @param img      Array 3d: all channels of the original image
#' @param p_mat    Matrix: the correction matrix of scaling factors
#' @return      Array 3d: the corrected image
correct_image <- function(img, p_mat) {
  img_new <- array(0, dim=dim(img))
  n_ch <- dim(img)[3]
  for (i in 1:n_ch) {
    for (j in 1:n_ch) {
      img_new[,,i] <- img_new[,,i] + img[,,j]*p_mat[j,i]
    }
  }
  img_new <- round(img_new)
  img_new[img_new < 0] <- 0
  return(img_new)
}