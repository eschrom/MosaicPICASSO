% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_functions.R
\name{correct_image}
\alias{correct_image}
\title{Apply Correction Matrix to an Image}
\usage{
correct_image(img, p_mat)
}
\arguments{
\item{img}{Array 3d: all channels of the original image}

\item{p_mat}{Matrix: the correction matrix of scaling factors}
}
\value{
\if{html}{\out{<div class="sourceCode">}}\preformatted{ Array 3d: the corrected image
}\if{html}{\out{</div>}}
}
\description{
Apply the correction matrix of scaling factors to every pair of channels to correct an image.
}
\details{
This function corrects a channel using a full correction matrix. Each channel is set to its original self,
minus a scaled version of each of the other channels. With respect to the correction matrix,
the column is the focal channel to be corrected, and the rows are the scaled factors for each of the other
channels to be subtracted.
}
