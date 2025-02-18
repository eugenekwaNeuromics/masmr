\name{checkReadyLoop}
\alias{checkReadyLoop}
\title{checkReadyLoop}
\description{
Check if ready for looping through images.
}
\usage{
checkReadyLoop()
}
\arguments{
}
\value{
  Boolean indicating if all checks have been passed. Warning messages returned if a test has failed.
}
\details{
}
\seealso{
}
\examples{
isReady <- checkReadyLoop()
if(!isReady){
  stop('Not ready to run loop!')
}
}
