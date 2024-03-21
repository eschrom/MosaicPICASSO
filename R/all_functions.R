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
  print("Loading image.")                                                       # Load img to 3d array format:
  IMG <- load_image(file)                                                       # x, y, c where c is channels
  N_CH <- dim(IMG)[3]                                                           # Number of channels
  if (enhance) {                                                                # If image will be contrast-enhanced
    frc_sat <- rep(NA, N_CH)                                                    # Record fraction of saturated pix
    for (i in 1:N_CH) {                                                         # in each channel (brightness = 255)
      frc_sat[i] <- sum(IMG[,,i] == 255) / prod(dim(IMG)[1:2])
    }
  }
  if (!is.null(groups)) {                                                       # If user specified groups
    print("Splitting image into channel groups for correction.")                # of channels to correct separately
    if (any(table(unlist(groups)) > 1)) {                                       # Make sure no channels are included
      stop("Error: Each channel can only be used in one group.")                # in multiple groups
    }                                                                           # Check to see if any channels are
    if (length(unlist(groups)) != N_CH) {                                       # not used at all - if so, warn user
      print("Warning: Channels that are not included in a group will not be used for correction.")
    }
    img <- vector("list", length(groups))                                       # Initialize a list of images &
    mat <- vector("list", length(groups))                                       # correction matrices, 1 per group
    names(mat) <- paste("G",1:length(groups))                                   # Give names to correction matrices
    for (g in 1:length(groups)) {                                               # Split image into several sub-
      img[[g]] <- IMG[,,groups[[g]],drop=F]                                     # images, each with a subset of
    }                                                                           # channels, stored in one list
  } else {                                                                      # If the user has not specified                                                               
    img <- vector("list", 1)                                                    # groups, still put the image and
    img[[1]] <- IMG                                                             # correction matrix in a list, but
    mat <- vector("list", 1)                                                    # only one element for each list
  }
  for (g in 1:length(img)) {                                                    # Correct each sub-image separately
    print(paste("Beginning channel group ", g, ".", sep=""))
    n_ch <- dim(img[[g]])[3]                                                    # Number of channels in this group
    p_mat_0 <- diag(n_ch)                                                       # Initial correction matrix                                                        
    p_mat_1 <- diag(n_ch)                                                       # Updated correction matrix
    print("Calculating correction matrix.")
    for (i in 1:cycles) {                                                       # For each cycle of correction
      print(paste("   Cycle ", i,".",sep=""))
      if (i > 1) {                                                              # If this is not the 1st cycle
        tmp_img <- correct_image(img[[g]], p_mat_0)                             # correct sub-image using past cycle
      } else {                                                                  # If this is the 1st cycle
        tmp_img <- img[[g]]                                                     # just start with original sub-image
      }                                                                         # Calculate new correction matrix
      tmp_mat <- calc_Pmat(tmp_img, method, q_thr, tile_size, stride, pix_ss, cores)# from current sub-image version
      p_mat_1 <- update_Pmat(p_mat_0, tmp_mat)                                  # Update correction matrix
      p_mat_0 <- p_mat_1                                                        # which then becomes new initial
    }
    print("Applying correction matrix.")
    img[[g]] <- correct_image(img[[g]], p_mat_1)                                # Correct the sub-image and store
    mat[[g]] <- p_mat_1                                                         # the final correction matrix.
  }
  if (!is.null(groups)) {                                                       # If user specified channel groups
    for (i in 1:N_CH) {                                                         # For each original channel
      grp <- which(unlist(lapply(groups, function(x,y) {y %in% x}, y=i)))       # Find which group its in
      idx <- which(groups[[grp]] == i)                                          # and which index w/in that group
      IMG[,,i] <- img[[grp]][,,idx]                                             # Reassemble full image with all
    }                                                                           # channels in correct order
  } else {                                                                      # If user didn't specify groups
    IMG <- img[[1]]                                                             # reassembly is unnecessary
  }
  if (enhance) {                                                                # If contrast-enhancing final image
    print("Enhancing contrast.")
    for (i in 1:N_CH) {                                                         # Enhance each channel separately
      if (frc_sat[i] == 0) {                                                    # If orig channel had no saturation
        mult <- 255                                                             # Use a multiplier of 255
      } else {                                                                  # Otherwise, multiplier is a quantile
        mult <- unname(quantile(c(IMG[,,i]), 1-frc_sat[i]))                     # of all pix values in new images
      }                                                                         # Multiplier =0 if new image is black
      if (mult > 0) {                                                           # If new image is not totally black
        IMG[,,i] <- round(IMG[,,i] * (255/mult))                                # Boost all pix by constant multiple
        IMG[,,i][IMG[,,i] > 255] <- 255                                         # Any pix >255 are set back to 255
      }
    }
  }
  print("Saving corrected image.")
  out_file <- paste(sub(".tif", "", file), "_NEW.tif", sep="")                  # Add _NEW to file name
  ijtiff::write_tif(IMG, out_file, overwrite=T)                                 # Save corrected image as .tif
  print("Saving correction matrix.")
  out_file <- paste(sub(".tif", "", file), "_MAT.xlsx", sep="")                 # Add _MAT to file name
  openxlsx::write.xlsx(mat, out_file)                                           # Save correction matrices, where
  print("Finished.")                                                            # matrix for each group is a sheet
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
  img_dim <- c(dim(img[[1]])[1], dim(img[[1]])[2], length(img))                 # Image is list of 2d channels
  for (i in 1:length(img)) {                                                    # For each channel
    img[[i]] <- img[[i]]*255                                                    # put pixel values on 0-255 scale
    mode(img[[i]]) <- "integer"                                                 # and store in integer format
  }
  img <- simplify2array(img)                                                    # Collapse channel list to 3d array:
  return(img)                                                                   # x, y, c where c is channels
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
  n_ch <- dim(img)[3]                                                           # Number of channels in image
  p_mat <- diag(n_ch)                                                           # Initialize correction matrix
  idx <- expand.grid(1:n_ch, 1:n_ch)                                            # All possible pairs of channels
  idx <- idx[idx[,1] != idx[,2], ]                                              # except if both channels are same
  img_chp <- chop_image(img, tile_size, stride)                                 # Chop image into tiles, listed by channel
  num_cores <- min(parallel::detectCores() - 1, cores)                          # Set up cluster for parallel 
  make_cluster <- parallel::makeCluster(num_cores, type = "PSOCK")              # for loop
  doParallel::registerDoParallel(cl = make_cluster)            
  `%dopar%` <- foreach::`%dopar%`                                               # For each pair of channels
  a <- foreach::foreach (i = c(1:nrow(idx)), .export=c('calc_scalar','calc_ssim','calc_MI','calc_corr')) %dopar% {
    calc_scalar(img_chp, idx[i,1], idx[i,2], method, q_thr, pix_ss)             # calculate correction scaling factor
  }                                                                             # returned as a list of numbers
  parallel::stopCluster(cl = make_cluster)                                      # Disassemble the parallel cluster
  a <- unlist(a)                                                                # Scaling factors as vector not list
  k <- 1                                                                        # Start with first scaling factor
  for (i in 1:nrow(idx)) {                                                      # For each pair of channels
    if (idx[i,1] != idx[i,2]) {                                                 # Make sure channels aren't same
      p_mat[idx[i,1], idx[i,2]] <- -a[k]                                        # Put scaling factor in its position
      k <- k + 1                                                                # in correction matrix, and move on
    }                                                                           # to the next scaling factor
  }
  return(p_mat)                                                                 # Return the full correction matrix
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
  coors <- expand.grid(row = seq(0,floor(dim(img)[1]/stride)-1)*stride + 1,     # Coordinates of top-left
                       col = seq(0,floor(dim(img)[2]/stride)-1)*stride + 1)     # corner of each tile
  blocks <- kronecker(matrix(1:nrow(coors), nrow=length(unique(coors$row)),     # (Nearly) full 2d extent of image
                             ncol=length(unique(coors$col)), byrow=F),          # with all pix labeled by which
                      matrix(1, stride, stride))                                # tile they belong to
  d1 <- dim(blocks)[1]                                                          # Width of block array
  d2 <- dim(blocks)[2]                                                          # Height of block array
  tiles <- vector("list", dim(img)[3])                                          # Initialize list of tile stacks
  for (i in 1:dim(img)[3]) {                                                    # For each channel
    tiles[[i]] <- lapply(split(img[1:d1,1:d2,i], blocks), matrix, nrow = stride)# Make a list of 2d tiles, then 
    tiles[[i]] <- simplify2array(tiles[[i]])                                    # convert to 3d array of tiles
    if (tile_size < stride) {                                                   # If tiles should be
      tiles[[i]] <- tiles[[i]][1:tile_size, 1:tile_size, ]                      # smaller, cut them
    }
  }                                                                             # Return list of 3d array of tiles
  return(tiles)                                                                 # xy = spatial position, z = tile ID
}                                                                               # List index = channel

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
  img_chp <- img_chp[c(i,j)]                                                    # 3d tile arrays for focal channels
  ssim <- calc_ssim(img_chp)                                                    # SSIM btw channels for each tile
  thr <- quantile(ssim, q_thr)                                                  # SSIM threshold to select tiles
  img_chp_sub <- vector("list", 2)                                              # Initialize a 2-element list
  for (i in 1:2) {                                                              # Fill in that list with only the
    img_chp_sub[[i]] <- img_chp[[i]][,,ssim <= thr]                             # tiles below threshold SSIM, for
  }                                                                             # both focal channels
  if (method == "MI") {                                                         # If minimization target is MI
    obj_func <- function(a, img1, img2) {                                       # define objective function to 
      return(calc_MI(img1, img2 - a*img1, pix_ss))                              # return MI btw ch1 & scaled ch2
    }
  } else if (method == "Cor") {                                                 # If minimization target is |cor|
    obj_func <- function(a, img1, img2) {                                       # define objective function to
      return(calc_corr(img1, img2 - a*img1, pix_ss))                            # return |cor| btw ch1 & scaled ch2
    }
  }
  a_try <- as.list(seq(0,1,by=0.01))                                            # Try a range of scaling factors &
  obj_out <- lapply(a_try, obj_func, img1=img_chp_sub[[1]], img2=img_chp_sub[[2]]) # list the MI or |cor| for each
  out <- unlist(a_try)[which.min(unlist(obj_out))]                              # Choose scaling factor value that
  if (length(out)==0) {                                                         # minimizes MI or |cor|
    return(0)                                                                   # If cor/MI is always NA, 
  } else {                                                                      # choose 0 as the scaling factor
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
  ss <- rep(NA, dim(tiles[[1]])[3])                                             # Consider each pair of tiles 
  for (i in 1:length(ss)) {                                                     # Measure SSIM between them
    ss[i] <- SpatialPack::SSIM(tiles[[1]][,,i], tiles[[2]][,,i])$SSIM
  }
  return(ss)                                                                    # Return vector of SSIM values for
}                                                                               # each tile pair

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
calc_corr <- function(img1, img2, pix_ss) {                                     # Melt 3d array of tiles to vector
  img1 <- as.vector(img1)                                                       # of pix values: values must match
  img2 <- as.vector(img2)                                                       # but spatial positions don't matter
  if (!is.character(pix_ss)) {                                                  # If a number of pix is given that's
    pix_ss <- min(length(img1), pix_ss)                                         # less than the number present 
    keep_pix <- sample(length(img1), pix_ss)                                    # randomly sample pixels & use only
    img1 <- img1[keep_pix]                                                      # these from each image
    img2 <- img2[keep_pix]
  }
  if (length(unique(img1)) == 0 || length(unique(img2))==0) {                   # If no variation in pix values
    return(NA)                                                                  # cannot compute correlation.
  } else {                                                                      # Otherwise, compute |cor| of pix
    cor_coef <- abs(cor(img1, img2))                                            # values across the two images
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
calc_MI <- function(img1, img2, pix_ss) {                                       # Melt 3d array of tiles to vector
  img1 <- as.vector(img1)                                                       # of pix values: values must match
  img2 <- as.vector(img2)                                                       # but spatial positions don't matter
  if (!is.character(pix_ss)) {                                                  # If a number of pix is given that's
    pix_ss <- min(length(img1), pix_ss)                                         # less than the number present 
    keep_pix <- sample(length(img1), pix_ss)                                    # randomly sample pixels & use only
    img1 <- img1[keep_pix]                                                      # these from each image
    img2 <- img2[keep_pix]
  }
  if (length(unique(img1)) == 0 || length(unique(img2))==0) {                   # If no variation in pix values
    return(NA)                                                                  # cannot compute correlation.
  } else {                                                                      # Otherwise, compute MI of pix
    img1 <- round(img1)                                                         # values across the two images
    img2 <- round(img2)                                                         # Make sure pix values are integers
    img2 <- img2 - min(img2)                                                    # No negative entries in image 2
    ent1 <- unname(table(img1))                                                 # Histogram of pix values for each
    ent2 <- unname(table(img2))                                                 # image
    ent1 <- ent1 / sum(ent1)                                                    # Transform histograms to
    ent2 <- ent2 / sum(ent2)                                                    # probabilities
    ent1 <- -sum(ent1 * log2(ent1))                                             # Calculate entropies from 
    ent2 <- -sum(ent2 * log2(ent2))                                             # probabilities
    imgb <- img1*(max(img2)+1) + img2                                           # Composite image with unique pix 
    entb <- unname(table(imgb))                                                 # values for every unique combo of
    entb <- entb / sum(entb)                                                    # pix values btw images 1 & 2
    entb <- -sum(entb * log2(entb))                                             # Histogram -> probs -> entropy
    mi <- ent1 + ent2 - entb                                                    # Calculate MI from entropies
    return((2*mi)/(ent1 + ent2))                                                # Return normalized MI
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
  return(p_mat_1 %*% p_mat_0)                                                   # Matrix multiplication to update
}                                                                               # correction matrix

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
  img_new <- array(0, dim=dim(img))                                             # Initialize new image as all black
  n_ch <- dim(img)[3]                                                           # Number of channels
  for (i in 1:n_ch) {                                                           # For each pair of channels,
    for (j in 1:n_ch) {                                                         # construct new image by adding 
      img_new[,,i] <- img_new[,,i] + img[,,j]*p_mat[j,i]                        # orig channel minus scaled 
    }                                                                           # versions of other channels
  }                                                                             # Correction matrix is 1 on diag
  img_new <- round(img_new)                                                     # and non-positive elsewhere
  img_new[img_new < 0] <- 0                                                     # Round pixel values & make sure
  return(img_new)                                                               # none are negative
}