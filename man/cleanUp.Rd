\name{cleanUp}
\alias{cleanUp}
\title{cleanUp}
\description{
Selective cleaning of environment.
}
\usage{
cleanUp( cleanSlate=c() )
}
\arguments{
\item{ cleanSlate }{ Character vector. Containing names of objects to save. }
}
\value{
  Cleaned global environment, while saving objects in provided list.
}
\examples{
establishCleanSlate()
cleanUp()

## Alternatively:
# cleanSlate <- establishCleanSlate()
# cleanUp(cleanSlate)
}
