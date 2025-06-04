\name{maxIntensityProject}
\alias{maxIntensityProject}
\title{maxIntensityProject}
\description{
Returns the maximum intensity projection of a Z-stack.
}
\usage{
maxIntensityProject(
  im,
  zDim = NULL
)
}
\arguments{
\item{ im }{ A 2-4D numeric array, or a list of 2-3D arrays. }
\item{ zDim }{ Integer. Specifies the dimension that Z planes are found in. When unspecified, defaults to expected axis (see Details). }
}
\value{
  A list of 2D matrices, or a single 2D matrix.
}
\details{
  Images are loaded with the following format assumed: [Number of Channels, Number of rows, Number of Columns, Number of Z slices]. See \code{readImage} for details.
  
  When a 4D array is provided, the first dimension is assumed to be Channel, and thus ignored. \code{zDim} is thus set to 3 as default.
  
  When provided as a list, it is assumed each element in the list covers a single channel and comprises a 3D array with the format [Number of rows, Number of Columns, Number of Z slices]. Again, \code{zDim} is thus set to 3 as default.
}
\seealso{
\link{readImage}
\link{readImageList}
}
\examples{
}
