\name{scoreClassifier}
\alias{scoreClassifier}
\title{scoreClassifier}
\description{
Returns score for a classifier.
}
\usage{
scoreClassifier(
    spotcalldf,
    labels,
    variablesExcluded,
    variablesIncluded,
    trainIndex,
    modelFormula,
    classifyFunction,
    LHS = 'OUT',
    params = get('params', envir = globalenv()),
    seed = 12345,
    predictType = 'response',
    returnModel = F,
    ...
)
}
\arguments{
\item{ spotcalldf }{ Dataframe containing per-pixel metrics. }
\item{ labels }{  A character indicating column in \code{spotcalldf} or a vector of row-wise identities. }
\item{ variablesExcluded }{  Character vector. Columns in \code{spotcalldf} to avoid in generating classifier model. If unspecified, automatically inferred. }
\item{ variablesIncluded }{  Character vector. Columns in \code{spotcalldf} to include for generating classifier model. If unspecified, automatically inferred. }
\item{ trainIndex }{ Boolean vector. Indicates which rows of \code{spotcalldf} are training data for classifier model. If unspecified, defaults to all rows. }
\item{ modelFormula }{ Character. Formula of format \code{ Y ~ X + ...}. If unspecified, inferred from \code{variablesIncluded} and \code{variablesExcluded}. }
\item{ classifyFunction }{ Function for classification. If unspecified, defaults to logisitic regression (see \code{glm}). }
\item{ LHS }{ Character. Left hand side of formula. }
\item{ seed }{ Integer. Seed for pseudo-random number generation. Defaults to \code{params$seed} if available. }
\item{ predictType }{ Character. For input into the \code{type} parameter of \code{predict} function: see \code{predict}. }
\item{ returnModel }{ Boolean. Determines if classifier model is returned. }
\item{ ... }{ Other parameters passed to \code{classifyFunction} }

}
\value{
  Numeric vector of scores per row in dataframe.
}
\details{
  If \code{returnModel} or if \code{predict} is incompatible with \code{classifyFunction}, the classifier model is returned (as \code{classifierModel} in global environment).

  If \code{classifyFunction} cannot take \code{modelFormula} as input, \code{modelFormula} is returned instead (as \code{modelFormula} in global environment).
}
\seealso{
\link{glm}
\link{predict}
}
\examples{
scores <- scoreClassifier(spotcalldf, 'BLANK')
spotcalldf$ProbabilityOfBeingBLANK <- scores

## Logistic regression is default, i.e.
## classifyFunction <- glm
## formals(classifyFunction)$family <- 'binomial'
## When predictType is 'response', probabilities are returned.

## Users may also wish to confirm model has converged (available for glm), e.g.:
scores <- scoreClassifier(spotcalldf, 'BLANK', returnModel = T) #classifierModel written to global environment
if( classifierModel$converged ){
  spotcalldf$ProbabilityOfBeingBLANK <- scores
}
}
